# -*- coding: utf-8 -*-

import os
import random
import math
import copy

import util
import pdbatoms
import v3
import data

# ---------------------------------------------------
# Standard units for this moduele:
# ---------------------------------------------------
# Velocity units: Ångstrom/ps
# Energy units: Da⋅Ångstrom^2/ps^2
# Mass units: Da
# ---------------------------------------------------

# 1 Da = 1.66E-27 kg
# boltzmann = 1.3806488E-23 J/K 
#           = 1.3806488E-23 kg⋅m^2/s^2/K 
#           = 1.3806488E-23/1.66E-27 Da⋅m^2/s^2/K
boltzmann_in_DaMSqPerSSqPerK = 8314.47148  

# momentum
# Da⋅Å/ps = 1.66E−27 kg⋅E-10 m/E-12 s = 1.66E-25 kg⋅m/s 
#         = 1.66E-25 N⋅s = 1.66E-25⋅E12⋅E12 pN⋅ps
momentum_DaAngPerPs_to_pNps = 1.66E-1

# energy/work units
# force = Da⋅Å/ps/ps = 1.66E-13 kg⋅m/s/s = 1.66E-13 N = 1.66E-1 pN
# work = energy = Da⋅Å⋅Å/ps/ps = Da⋅Å/ps/ps ⋅ Å = 1.66E-1 pNÅ
work_DaAngSqPerPsSq_to_pNAng = 1.66E-1
work_Nm_to_kcal = 0.000238846
avogadro = 6.02214179E23
work_pNAng_to_kcalPerMol = 1E-12*1E-10*work_Nm_to_kcal*avogadro

# velcoity units
vel_mPerS_to_AngsPerPs = 0.01 # from m/s to angstroms/ps
velsq_mPerS_to_AngsPerPs = 1.0E-4 # (m/s)^2 to (angstroms/ps)^2

# typical molecular-dynamics integration time-step
timestep_in_ps = 0.001 


# Atom maniuplation

def get_atoms_of_residues(soup, residues, selection=None):
  atoms = []
  for i in residues:
    res_atoms = soup.residue(i).atoms()
    if selection:
      res_atoms = [a for a in res_atoms if a.type in selection]
    atoms.extend(res_atoms)
  return atoms
  
def average_vel(atoms):
  momentum = v3.vector()
  mass = 0.0
  for a in atoms:
    momentum += v3.scale(a.vel, a.mass)
    mass += a.mass
  return v3.scale(momentum, 1.0/mass)

def add_vel_to_atoms(atoms, vel_diff):
  for a in atoms:
    a.vel_last = a.vel
    a.vel += vel_diff
    a.work_delta = \
        v3.dot(vel_diff, a.vel) * a.mass * timestep_in_ps * \
        work_DaAngSqPerPsSq_to_pNAng


# Heating functions

def maxwell_velocity(temp, mass):
  """
  input: temp in K, mass in Da
  output: velocity in Ångstroms/ps
  """
  velsq_ave = boltzmann_in_DaMSqPerSSqPerK * temp / mass # in (m/s)^2
  return random.gauss(0, math.sqrt(velsq_ave)) * vel_mPerS_to_AngsPerPs


def mean_energy(temp, n_degree_of_freedom):
  "output in Da (angs/ps)^2"
  return n_degree_of_freedom * 0.5 * boltzmann_in_DaMSqPerSSqPerK * \
         temp * velsq_mPerS_to_AngsPerPs # Da (angstroms/ps)^2


def random_energy(temp, n_degree_of_freedom):
  average = mean_energy(temp, n_degree_of_freedom);
  std_dev = math.sqrt(average)
  return random.gauss(average, std_dev)

  
def kinetic_energy(atoms):
  en = 0.0
  for a in atoms:
    vel = v3.mag(a.vel)
    en += 0.5 * a.mass * vel * vel
  return en


def gas_randomize(atoms, temp):
  "Heats residues uses gas approximation: vel: angstroms/ps"
  for atom in atoms:
    v3.set_vector(
        atom.vel,
        maxwell_velocity(temp, atom.mass),
        maxwell_velocity(temp, atom.mass),
        maxwell_velocity(temp, atom.mass))


def anderson_velocity_scale(atoms, temp, n_degree_of_freedom):
  "Scales velocity of atoms to energy at temp. Vel: angstroms/ps"
  target_energy = mean_energy(temp, n_degree_of_freedom)
  kin = kinetic_energy(atoms)
  if v3.is_similar_mag(kin, 0):
    gas_randomize(atoms, temp)
  else:
    scaling_factor = math.sqrt(target_energy / kin)
    for atom in atoms:
      v3.set_vector(atom.vel, v3.scale(atom.vel, scaling_factor))


# Push functinos

class PushApartByVel():
  """
  This is probably the most common case, push apart two domains
  using a maximum target velocity.
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
    self.soup = soup
    self.setup_domains()
    self.change_vels()
    self.calculate_output()
    self.append_output()


def make_puff_fn(
    domain1, domain2, target_val, dt=0.1, temp=None, 
    is_backbone_only=False, is_first_domain_only=False, 
    force_fname='md.puff.out'):
  strategy = PushApartByVel(
      domain1, domain2, target_val, dt, temp, 
      is_backbone_only, is_first_domain_only, 
      force_fname)
  pulse_fn = lambda soup: strategy.apply(soup)
  return pulse_fn


def read_puff_out(md_dir):
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
  The simplest push: a constant force, acceleartion taken from target_val.
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
  strategy = PushApartByAcc(
      domain1, domain2, target_val, dt, temp, 
      is_backbone_only, is_first_domain_only, force_fname)
  pulse_fn = lambda soup: strategy.apply(soup)
  return pulse_fn


# ATD sidechain heating

def make_atd_fn(i_residue, heating_temp, backbone_atoms):
  def gas_heat_sidechain(
      soup, i_residue, heating_temp, backbone_atoms):
    atoms = [a for a in soup.residue(i_residue).atoms() 
             if a.type not in backbone_atoms]
    gas_randomize(atoms, heating_temp)
  return lambda soup: gas_heat_sidechain(
      soup, i_residue, heating_temp, backbone_atoms)


# General rotation functions

# Rotational Units
# --------------------------------------------------------
# rotational-velocity = º/ps = E+12 º/s 
# rotational-acceleration = º/ps/ps = E+24º/s/s
# moment-of-inertia = Da⋅Å⋅Å 
#                   = 1.66E-27 kg⋅E-10m⋅E-10m 
#                   = 1.66E-47 kg⋅m^2
# torque = moment-of-inertia*rotational-acceleration
#        = Da⋅Å⋅Å ⋅ º/ps/ps 
#        = 1.66E-47 kg⋅m^2 ⋅ E+24⋅º/s/s
#        = 1.66E-23 º⋅m⋅kg⋅m/s/s 
#        = 1.66E-23 m⋅N
# force = torque/radius
#       = Da⋅Å⋅Å⋅º/ps/ps / Å 
#       = 1.66E-23 m⋅N/E-10m
#       = 1.66E-13 N
#       = 1.66 E-1 pN

def moment_of_inertia(atom, axis, anchor):
  "Returns moment in Da*angstroms^^2"
  r = atom.pos - anchor
  r_perp = v3.perpendicular(r, axis)
  r_len = v3.mag(r_perp)
  return atom.mass * r_len * r_len


def rotational_velocity(atom, axis, anchor):
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
  

def total_moment_of_inertia(atoms, axis, anchor):
  moments = [moment_of_inertia(atom, axis, anchor)
             for atom in atoms]
  return sum(moments)
    

def weighted_rotational_velocity(atoms, axis, anchor):
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
  

# Sidechain rotation functions

chi_topology = {
  'ARG': [ ['N', 'CA', 'CB', 'CG'],
           ['CA', 'CB', 'CG', 'CD'],
           ['CB', 'CG', 'CD', 'NE'],
           ['CG', 'CD', 'NE', 'CZ']],
  'ASN': [['N', 'CA', 'CB', 'CG'], ['CA', 'CB', 'CG', 'OD1']],
  'ASP': [['N', 'CA', 'CB', 'CG'], ['CA', 'CB', 'CG', 'OD1']],
  'CYS': [['N', 'CA', 'CB', 'SG']],
  'GLN': [ ['N', 'CA', 'CB', 'CG'],
           ['CA', 'CB', 'CG', 'CD'],
           ['CB', 'CG', 'CD', 'OE1']],
  'GLU': [ ['N', 'CA', 'CB', 'CG'],
           ['CA', 'CB', 'CG', 'CD'],
           ['CB', 'CG', 'CD', 'OE1']],
  'HIS': [['N', 'CA', 'CB', 'CG'], ['CA', 'CB', 'CG', 'ND1']],
  'ILE': [['N', 'CA', 'CB', 'CG1'], ['CA', 'CB', 'CG1', 'CD1']],
  'LEU': [['N', 'CA', 'CB', 'CG'], ['CA', 'CB', 'CG', 'CD1']],
  'LYN': [ ['N', 'CA', 'CB', 'CG'],
           ['CA', 'CB', 'CG', 'CD'],
           ['CB', 'CG', 'CD', 'CE'],
           ['CG', 'CD', 'CE', 'NZ']],
  'LYP': [ ['N', 'CA', 'CB', 'CG'],
           ['CA', 'CB', 'CG', 'CD'],
           ['CB', 'CG', 'CD', 'CE'],
           ['CG', 'CD', 'CE', 'NZ']],
  'LYS': [ ['N', 'CA', 'CB', 'CG'],
           ['CA', 'CB', 'CG', 'CD'],
           ['CB', 'CG', 'CD', 'CE'],
           ['CG', 'CD', 'CE', 'NZ']],
  'MET': [ ['N', 'CA', 'CB', 'CG'],
           ['CA', 'CB', 'CG', 'SD'],
           ['CB', 'CG', 'SD', 'CE']],
  'PHD': [['N', 'CA', 'CB', 'CG'], ['CA', 'CB', 'CG', 'OD1']],
  'PHE': [['N', 'CA', 'CB', 'CG'], ['CA', 'CB', 'CG', 'CD1']],
  'PRO': [ ['N', 'CA', 'CB', 'CG'],
           ['CA', 'CB', 'CG', 'CD'],
           ['CB', 'CG', 'CD', 'N'],
           ['CG', 'CD', 'N', 'CA']],
  'SER': [['N', 'CA', 'CB', 'OG']],
  'THR': [['N', 'CA', 'CB', 'OG1']],
  'TRP': [['N', 'CA', 'CB', 'CG'], ['CA', 'CB', 'CG', 'CD1']],
  'TYR': [['N', 'CA', 'CB', 'CG'], ['CA', 'CB', 'CG', 'CD1']],
  'VAL': [['N', 'CA', 'CB', 'CG1']]}


def get_res_chi_topology(res):
  """
  Gets the atom types for each chi angle in a residue type
  """
  res_type = res.type
  if res_type not in chi_topology:
    return []
  if res_type in ["HSE"]:
    res_type = "HIS"
  result = copy.deepcopy(chi_topology[res_type])
  if res_type == "ILE" and res.has_atom("CD1"):
    for i in range(0, len(result)):
      for j in range(0, len(result[i])):
        if result[i][j] == "CD":
          result[i][j] = "CD1"
          break
  return result


def get_n_chi(res):
  if chi_topology.has_key(res.type):
    return len(chi_topology[res.type])
  return 0


def calculate_chi(res, j):
  res_chi_topology = get_res_chi_topology(res)
  if j < len(res_chi_topology):
    p = [res.atom(atom_type).pos for atom_type in res_chi_topology[j]]
    return v3.normalize_angle(v3.dihedral(p[0], p[1], p[2], p[3]))
  raise ValueError, "No Chi%d angle for res %d" % (j, i)


def get_axis_anchor(res, i_chi):
  chi_topology = get_res_chi_topology(res)
  p = [res.atom(a).pos for a in chi_topology[i_chi]]
  axis = p[2] - p[1]
  anchor = p[2]
  return axis, anchor
  
    
def atoms_affected_by_chi(res, i_chi):
  """
  Identifies chi angle from atom_type: -1=backbone, 0=chi1, etc.
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
  return [a for a in res.atoms() if sidechain_nesting(a.type) >= i_chi]


def get_rot_vel_chi(res, i):
  axis, anchor = get_axis_anchor(res, i)    
  atoms = atoms_affected_by_chi(res, i)
  return weighted_rotational_velocity(atoms, axis, anchor)


def get_random_chi_rot_vel(res, i, temp):
  axis, anchor = get_axis_anchor(res, i)
  atoms = atoms_affected_by_chi(res, i)
  moment = total_moment_of_inertia(atoms, axis, anchor)
  energy = random_energy(temp, 3*len(atoms))
  return math.sqrt(2 * energy / moment)


def add_rot_vel_to_chi(res, i, target_rot_vel):
  axis, anchor = get_axis_anchor(res, i)
  atoms = atoms_affected_by_chi(res, i)
  add_rotational_velocity(atoms, target_rot_vel, axis, anchor)


class Rip:
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
  rip = Rip(i_res, heating_temp)
  return lambda soup: rip.apply(soup)

