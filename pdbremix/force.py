# encoding: utf-8

__doc__ = """
This module provides functions to calculate velocity changes
in Soup objects due to different types of applied forces. This
can be saved to restart files for AMBER, NAMD and GROMACS.

There are four functions that are to be used:

1. make_atd_fn(i_residue, heating_temp, backbone_atoms)

2. make_puff_fn(
    domain1, domain2, target_val, dt=0.1, temp=None, 
    is_backbone_only=False, is_first_domain_only=False, 
    force_fname='md.puff.out')

3. make_puff_acc_fn(
    domain1, domain2, target_val, dt=0.1, temp=None, 
    is_backbone_only=False, is_first_domain_only=False, 
    force_fname='md.puff.out')

4. make_rip_fn(i_res, heating_temp)


These are all function factories to generate a function in the
form:

   def pulse_fn(soup): 
     # change soup velocities
     
The calculations of velocities and energies used here are 
mainly:

- velocity: Ångstrom/ps
- work/energy: Da*angs/ps^2
- mass: Da

However, various conversions are needed to connect with various
thermostatistical identites as well as for output. For instance,
kinetic energies calculated from mass and velocities result in
the units: Da⋅Ångstrom^2/ps^2, which need to be converted to
kcal/mol for output. 
"""


import os
import random
import math
import copy

import util
import pdbatoms
import v3
import data



##########################################################
# Heating functions



def average_vel(atoms):
  """
  Returns the mass-averaged velocity of atoms.
  """
  momentum = v3.vector()
  mass = 0.0
  for a in atoms:
    momentum += v3.scale(a.vel, a.mass)
    mass += a.mass
  return v3.scale(momentum, 1.0/mass)


# Work/energy conversions
# force = Da⋅Å/ps/ps 
#       = 1.66E-13 kg⋅m/s/s 
#       = 1.66E-13 N 
#       = 1.66E-1 pN
# work/energy = Da⋅Å⋅Å/ps/ps 
#             = Da⋅Å/ps/ps ⋅ Å 
#             = 1.66E-1 pNÅ
work_DaAngSqPerPsSq_to_pNAng = 1.66E-1
work_Nm_to_kcal = 0.000238846
avogadro = 6.02214179E23
work_pNAng_to_kcalPerMol = 1E-12*1E-10*work_Nm_to_kcal*avogadro

# molecular-dynamics integration time-step
timestep_in_ps = 0.001 

def add_vel_to_atoms(atoms, vel_diff):
  """
  Adds vel_diff to the vel vector of atoms.
  """
  for a in atoms:
    a.vel_last = a.vel
    a.vel += vel_diff
    a.work_delta = \
        v3.dot(vel_diff, a.vel) * a.mass * timestep_in_ps * \
        work_DaAngSqPerPsSq_to_pNAng


# Velocity conversions
vel_mPerS_to_AngsPerPs = 0.01 
velsq_mPerS_to_AngsPerPs = 1.0E-4 
# Boltzmann constant
#    = 1.3806488E-23 J/K 
#    = 1.3806488E-23 kg ⋅ m^2/s^2/K 
#    = 1.3806488E-23 1/1.66E-27Da ⋅ m^2/s^2/K
#    = 8314.47148 Da⋅m^2/s^2/K
boltzmann_in_DaMSqPerSSqPerK = 8314.47148  

def maxwell_velocity(temp, mass):
  """
  Returns a velocity (in angs/ps) sampled from a Maxwell velocity
  distribution determined by mass and temp. 
  """
  velsq_ave = boltzmann_in_DaMSqPerSSqPerK * temp / mass 
  return random.gauss(0, math.sqrt(velsq_ave)) * \
         vel_mPerS_to_AngsPerPs


def mean_energy(temp, n_degree_of_freedom):
  """
  Returns the average energy (Da*angs/ps^2) of n degree of 
  freedom at temperature. if n_degree_of_freedom = 3,
  this is the average energy of a point particle.
  """
  return 0.5 * n_degree_of_freedom * temp * \
         boltzmann_in_DaMSqPerSSqPerK * \
         velsq_mPerS_to_AngsPerPs 


def random_energy(temp, n_degree_of_freedom):
  """
  Returns an energy (Da*angs/ps^2) sampled from a Maxwellian
  distribution of energies at temp.
  """
  average = mean_energy(temp, n_degree_of_freedom);
  std_dev = math.sqrt(average)
  return random.gauss(average, std_dev)

  
def kinetic_energy(atoms):
  """
  Returns the kinetic energy (Da*angs/ps^2) of the atoms.
  """
  en = 0.0
  for a in atoms:
    vel = v3.mag(a.vel)
    en += 0.5 * a.mass * vel * vel
  return en


def gas_randomize(atoms, temp):
  """
  Randomly assigns a velocity to atoms based on a Maxwellian
  distribution at temp.
  """
  for atom in atoms:
    v3.set_vector(
        atom.vel,
        maxwell_velocity(temp, atom.mass),
        maxwell_velocity(temp, atom.mass),
        maxwell_velocity(temp, atom.mass))


def make_atd_fn(i_residue, heating_temp, backbone_atoms):
  """
  Returns a function that applies velocity changes for the
  original Anisotropic Thermal Diffusion technique, which
  is essentially a localized thermostat on a single sidechain.
  """
  def gas_heat_sidechain(
      soup, i_residue, heating_temp, backbone_atoms):
    atoms = [a for a in soup.residue(i_residue).atoms() 
             if a.type not in backbone_atoms]
    gas_randomize(atoms, heating_temp)
  return lambda soup: gas_heat_sidechain(
      soup, i_residue, heating_temp, backbone_atoms)


def anderson_velocity_scale(atoms, temp, n_degree_of_freedom):
  """
  Classic Anderson thermometer: controls the temperature
  by scaling the velocity of atoms until the average energy
  reaches the average velocity consistent with the given temp.
  """
  target_energy = mean_energy(temp, n_degree_of_freedom)
  kin = kinetic_energy(atoms)
  if v3.is_similar_mag(kin, 0):
    gas_randomize(atoms, temp)
  else:
    scaling_factor = math.sqrt(target_energy / kin)
    for atom in atoms:
      v3.set_vector(atom.vel, v3.scale(atom.vel, scaling_factor))



##########################################################
# Pushing functions


def get_atoms_of_residues(soup, res_indices, atom_types=None):
  """
  Return atoms of soup that belong to residues indicated by
  res_indices and in atom_types.
  """
  atoms = []
  for i in res_indices:
    res_atoms = soup.residue(i).atoms()
    if atom_types:
      res_atoms = [a for a in res_atoms if a.type in atom_types]
    atoms.extend(res_atoms)
  return atoms
  

class PushApartByVel():
  """
  Strategy to appy forces to push apart two domains in a Soup.

  This class will not be directly used. Rather the associated
  make_puff_fn() will initialize an instance of PushApartByVel in
  a closure and return the apply(soup) method via a lambda
  to the pulse simulation of simulate.py. This closure will
  allow complex initiations to be remembered during the course
  of a pulsed simulation.
  """
  def __init__(
      self, domain1, domain2, target_val, dt=0.1,
      temp=None, is_backbone_only=False, 
      is_first_domain_only=True,
      force_fname='md.puff.out'):
    self.domain1 = domain1
    self.domain2 = domain2
    self.dt = dt
    self.target_val = target_val
    self.temp = temp
    self.force_fname = os.path.abspath(force_fname)
    self.is_backbone_only = is_backbone_only  
    self.is_first_domain_only = is_first_domain_only  

  def setup_domains(self):
    # scale the temperature (necessary at high pulling speeds)
    if self.temp:
      atoms = self.soup.atoms()
      anderson_velocity_scale(atoms, self.temp, 3*len(atoms))

    # select the atoms from the domains definition
    selection = data.backbone_atoms if self.is_backbone_only else None
    self.atoms1 = get_atoms_of_residues(self.soup, self.domain1, selection)
    self.atoms2 = get_atoms_of_residues(self.soup, self.domain2, selection)

    # get direction vectors based on domains
    self.disp2to1 = pdbatoms.get_center(self.atoms1) \
                  - pdbatoms.get_center(self.atoms2)
    self.axis2to1 = v3.norm(self.disp2to1)

    # calculate relative velocities between domains
    self.vel2to1 = average_vel(self.atoms1) - average_vel(self.atoms2)
    self.axis_vel2to1 = v3.parallel(self.vel2to1, self.axis2to1)
    self.vel = v3.dot(self.axis_vel2to1, self.axis2to1)

    if self.is_first_domain_only:
      self.move_atoms = self.atoms1
    else:
      self.move_atoms = self.atoms1 + self.atoms2

  def change_vels(self):
    # calculate the vel diff vector 
    target_axis_vel2to1 = v3.scale(self.axis2to1, self.target_val)
    diff_axis_vel2to1 = target_axis_vel2to1 - self.axis_vel2to1
    self.vel_diff = v3.dot(diff_axis_vel2to1, self.axis2to1)

    # now change velocities of movable atoms
    self.old_kinetic_energy = kinetic_energy(self.move_atoms)
    if self.is_first_domain_only:
      add_vel_to_atoms(self.atoms1, diff_axis_vel2to1)
    else:
      # apply half of vel_diff to each domain
      v3.scale(diff_axis_vel2to1, 0.5)
      add_vel_to_atoms(self.atoms1, diff_axis_vel2to1)
      add_vel_to_atoms(self.atoms2, -diff_axis_vel2to1)
    self.kinetic_energy = kinetic_energy(self.move_atoms)

  def calculate_output(self):
    work_applied = self.kinetic_energy - self.old_kinetic_energy
    work_applied *= work_DaAngSqPerPsSq_to_pNAng
    work_applied *= work_pNAng_to_kcalPerMol
    self.output_dict = {
      'separation': v3.mag(self.disp2to1),
      'mass': sum(a.mass for a in self.move_atoms),
      'target_vel': self.target_val,
      'vel': self.vel,
      'vel_diff': self.vel_diff,
      'work_applied': work_applied
    }

  def append_output(self):
    with open(self.force_fname, 'a') as f:
      out_s = str(self.output_dict)
      if not out_s.endswith('\n'):
        out_s += '\n'
      f.write(out_s)

  def apply(self, soup):
    """
    Apply velocity changes to the soup that will induce the
    pulling. This is the key method that will be exported to the
    pulse() function of simulate.py.
    """
    self.soup = soup
    self.setup_domains()
    self.change_vels()
    self.calculate_output()
    self.append_output()


def make_puff_fn(
    domain1, domain2, target_val, dt=0.1, temp=None, 
    is_backbone_only=False, is_first_domain_only=False, 
    force_fname='md.puff.out'):
  """
  Returns a function for simulate.py that implements 
  PushApartByVel. 

  This is the main interface to the pulsing function.
  """
  strategy = PushApartByVel(
      domain1, domain2, target_val, dt, temp, 
      is_backbone_only, is_first_domain_only, 
      force_fname)
  pulse_fn = lambda soup: strategy.apply(soup)
  return pulse_fn


def read_puff_out(md_dir):
  """
  Yields a dictionary representing the properties of 
  PushApartByVel at each frame of a pulsed simulation from the
  specified md.puff.out file.
  """
  # get time in ps, typical MD step is 0.001 ps = 1 fs
  config = os.path.join(md_dir, 'md.puff.config')
  parms = util.read_dict(config)
  dt = 0.001*parms['n_step_per_pulse']
  time = 0.0
  for line in open(os.path.join(md_dir, 'md.puff.out')):
    entry = eval(line)
    entry['time'] = time
    yield entry
    time += dt


class PushApartByAcc(PushApartByVel):
  """
  A Strategy class that push apart two domains in Soup
  with constant acceleartion.
  """
  def change_vels(self):
    # calculate the vel diff vector to apply
    diff_vel = self.target_val * self.dt
    diff_axis_vel2to1 = v3.scale(self.axis2to1, diff_vel)
    self.vel_diff = v3.dot(diff_axis_vel2to1, self.axis2to1)
    if self.is_first_domain_only:
      add_vel_to_atoms(self.atoms1, diff_axis_vel2to1)
    else:
      # apply half of vel_diff to each domain
      diff_axis_vel2to1 = v3.scale(diff_axis_vel2to1, 0.5)
      add_vel_to_atoms(self.atoms1, diff_axis_vel2to1)
      add_vel_to_atoms(self.atoms2, -diff_axis_vel2to1)

  def calculate_output(self):
    PushApartByVel.calculate_output(self)
    self.output_dict['target_acc'] = self.target_val
    del self.output_dict['target_vel']


def make_puff_acc_fn(
    domain1, domain2, target_val, dt=0.1, temp=None, 
    is_backbone_only=False, is_first_domain_only=False, 
    force_fname='md.puff.out'):
  """
  Returns a function that implements PushApartByAcc.
  """
  strategy = PushApartByAcc(
      domain1, domain2, target_val, dt, temp, 
      is_backbone_only, is_first_domain_only, force_fname)
  pulse_fn = lambda soup: strategy.apply(soup)
  return pulse_fn


##########################################################
# 3. Rotation functions

# Rotational Units
# --------------------------------------------------------
# rotational-velocity:
#    1 º/ps = E+12 º/s 
# rotational-acceleration:
#    1 º/ps/ps = E+24º/s/s
# moment-of-inertia
#    1 Da⋅Å⋅Å = 1.66E-27 kg⋅E-10m⋅E-10m 
#             = 1.66E-47 kg⋅m^2
# torque = moment-of-inertia*rotational-acceleration
#    Da⋅Å⋅Å⋅º/ps/ps  = 1.66E-47 kg⋅m^2 ⋅ E+24⋅º/s/s
#                   = 1.66E-23 º⋅m⋅kg⋅m/s/s 
#                   = 1.66E-23 º⋅m⋅N
# force = torque/radius
#    Da⋅Å⋅Å⋅º/ps/ps / Å = 1.66E-23 m⋅N/E-10m
#                      = 1.66E-13 N
#                      = 1.66 E-1 pN


def moment_of_inertia(atom, axis, anchor):
  """
  Returns the moment (DaAng^^2) of the atom connected to
  anchor around axis.
  """
  r = atom.pos - anchor
  r_perp = v3.perpendicular(r, axis)
  r_len = v3.mag(r_perp)
  return atom.mass * r_len * r_len


def total_moment_of_inertia(atoms, axis, anchor):
  """
  Returns the total moment (DaAng^^2) of a bunch of atoms that
  are connected to anchor around axis.
  """
  moments = [moment_of_inertia(atom, axis, anchor)
             for atom in atoms]
  return sum(moments)
    

def rotational_velocity(atom, axis, anchor):
  """
  Returns the rotational velocity (rad/ps) of the atom connected
  to anchor around axis.
  """
  r = atom.pos - anchor
  r_perp = v3.perpendicular(r, axis)
  vel_perp = v3.perpendicular(atom.vel, axis)
  vel_tang = v3.perpendicular(vel_perp, r_perp)
  pos_ref = v3.cross(axis, r_perp)
  if v3.dot(vel_tang, pos_ref) < 0.0:
    sign = -1.0
  else:
    sign = 1.0
  if v3.is_similar_mag(v3.mag(r_perp), 0):
    result = 0.0
  else:
    result = sign * v3.mag(vel_tang) / v3.mag(r_perp)
  return result
  

def weighted_rotational_velocity(atoms, axis, anchor):
  """
  Returns the average rotational velocity (rad/ps) of a bunch of
  a atoms weighted by the moment of the atom.
  """
  moments = [ \
      moment_of_inertia(atom, axis, anchor) for atom in atoms]
  total_moment = sum(moments)
  weights = [moment / total_moment for moment in moments]
  rot_vels = [ \
      rotational_velocity(atom, axis, anchor) for atom in atoms]
  weighted_rot_vels = [ \
      rot_vel*weight for rot_vel, weight in zip(rot_vels, weights)]
  return sum(weighted_rot_vels)


def add_rotational_velocity(atoms, rot_vel, axis, anchor):
  """
  Adds the rot_vel to the vel vector of atoms with respect
  to the rotation around axis and attached to anchor.
  """
  for atom in atoms:
    r_perp = v3.perpendicular(atom.pos - anchor, axis)
    v_tang_dir = v3.cross(axis, r_perp)
    v_tang_dir_len = v3.mag(v_tang_dir)
    if v3.is_similar_mag(v_tang_dir_len, 0):
      v_tang = v3.vector()
    else:
      v_new_len = rot_vel * v3.mag(r_perp)
      v_tang = v3.scale(v_tang_dir, v_new_len/v_tang_dir_len)
    atom.vel += v_tang
  

##########################################################
# Sidechain rotation functions

def get_res_chi_topology(res):
  """
  Returns the chi topology for a given res, which is a list of
  atoms that are affected if one rotates the chi0, chi1... 
  dihedral angle.
  """
  res_type = res.type
  if res_type not in data.chi_topology:
    return []
  if res_type in ["HSE"]:
    res_type = "HIS"
  result = copy.deepcopy(data.chi_topology[res_type])
  if res_type == "ILE" and res.has_atom("CD1"):
    for i in range(0, len(result)):
      for j in range(0, len(result[i])):
        if result[i][j] == "CD":
          result[i][j] = "CD1"
          break
  return result


def get_n_chi(res):
  """
  Returns the number of chi angles that res has.
  """
  if data.chi_topology.has_key(res.type):
    return len(data.chi_topology[res.type])
  return 0


def calculate_chi(res, j):
  """
  Returns the angle for the j'th chi dihedral angle of res.
  """
  res_chi_topology = get_res_chi_topology(res)
  if j < len(res_chi_topology):
    p = [res.atom(atom_type).pos for atom_type in res_chi_topology[j]]
    return v3.normalize_angle(v3.dihedral(p[0], p[1], p[2], p[3]))
  raise ValueError, "No Chi%d angle for res %d" % (j, i)


def get_axis_anchor(res, i_chi):
  """
  Returns the axis of rotation and an anchor point of the i_chi
  chi dihedral angle of res.
  """
  chi_topology = get_res_chi_topology(res)
  p = [res.atom(a).pos for a in chi_topology[i_chi]]
  axis = p[2] - p[1]
  anchor = p[2]
  return axis, anchor
  
    
def atoms_affected_by_chi(res, i_chi):
  """
  Returns the atoms in res that will be rotated if the i_chi
  chi dihedral angle is rotated.
  """
  def sidechain_nesting(atom_type):
    label = atom_type
    while label[-1].isdigit():
      label = label[:-1]
    while label[0].isdigit():
      label = label[1:]
    if len(label) < 2:
      nesting = -1
    else:
      nesting = "ABGDEZH".find(label[1]) - 2
      if label[0] == "H":
        nesting += 1
    return nesting
  return [a for a in res.atoms()
          if sidechain_nesting(a.type) >= i_chi]


def get_rot_vel_chi(res, i_chi):
  """
  Returns the weighted rotational velocity of the atoms 
  that are rotated by the ith chi angle.
  """
  axis, anchor = get_axis_anchor(res, i_chi)    
  atoms = atoms_affected_by_chi(res, i_chi)
  return weighted_rotational_velocity(atoms, axis, anchor)


def get_random_chi_rot_vel(res, i, temp):
  """
  Returns a random energy from a Maxwellian energy distribution
  consistent with the number of degrees of freedom of the
  atoms involved in the i_chi angle.
  """
  axis, anchor = get_axis_anchor(res, i)
  atoms = atoms_affected_by_chi(res, i)
  moment = total_moment_of_inertia(atoms, axis, anchor)
  energy = random_energy(temp, 3*len(atoms))
  return math.sqrt(2 * energy / moment)


def add_rot_vel_to_chi(res, i, target_rot_vel):
  """
  Adds target_rot_vel to the atoms affected by i around the
  chi axis.
  """
  axis, anchor = get_axis_anchor(res, i)
  atoms = atoms_affected_by_chi(res, i)
  add_rotational_velocity(atoms, target_rot_vel, axis, anchor)


class Rip:
  """
  Strategy class to implement the rotational velocity changes
  to a soup object. To be used by make_rip_fn.
  """
  def __init__(self, i_res, heating_temp):
    self.i_res = i_res
    self.heating_temp = heating_temp
    self.mean_chis = None
    self.max_delta_chi = v3.radians(60)

  def apply(self, soup):
    res = soup.residue(self.i_res)
    atoms = res.atoms()
    n_chi = get_n_chi(res)

    if self.mean_chis is None:
      self.mean_chis = [calculate_chi(res, i) for i in range(n_chi)]

    rot_vels = [get_rot_vel_chi(res, i) for i in range(n_chi)]
    for atom in atoms:
      v3.set_vector(atom.vel, 0.0, 0.0, 0.0)

    for i_chi in reversed(range(n_chi)):
      chi = calculate_chi(res, i_chi)
      delta_chi = v3.normalize_angle(chi - self.mean_chis[i])
      target_rot_vel = get_random_chi_rot_vel(
          res, i_chi, self.heating_temp)
      if abs(delta_chi) > self.max_delta_chi:
        if delta_chi > self.max_delta_chi:
          target_rot_vel = -target_rot_vel
      else:
        if rot_vels[i_chi] < 0.0:
          target_rot_vel *= -target_rot_vel
      add_rot_vel_to_chi(res, i_chi, target_rot_vel)

    anderson_velocity_scale(atoms, self.heating_temp, 3*len(atoms))


def make_rip_fn(i_res, heating_temp):
  """
  Returns a function for simulate.py that implements 
  Rip. 
  """
  rip = Rip(i_res, heating_temp)
  return lambda soup: rip.apply(soup)

