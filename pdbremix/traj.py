import os
import copy

import util
import v3
import pdbatoms
import data
import namd
import amber
import gromacs


def open_trajectory(name):
  if os.path.isfile(name + '.trr'):
    print "Detected GROMACS trajectory"
    return gromacs.Trajectory(name)
  if os.path.isfile(name + '.psf'):
    print "Detected NAMD trajectory"
    return namd.Trajectory(name)
  if os.path.isfile(name + '.top'):
    print "Detected AMBER trajectory"
    return amber.Trajectory(name)
  raise ValueError("No trajectories found at all!")


def make_equilibrium_pdb(trj, out_name):
  """
  Calculates the equilibrium structure from the
  last half of a given trajectory trj
  """
  start = trj.n_frame // 2
  end = trj.n_frame
  
  if start < 0 or start >= trj.n_frame:
    raise IndexError("Start frame out of range")
  if end < 0 or end > trj.n_frame+1:
    raise IndexError("End frame out of range")
  
  sum_soup = trj.soup.copy()
  
  for a in sum_soup.atoms():
    a.pos.set(0, 0, 0)
  
  n_frame = 0
  for i in range(start, end):
    trj.load_frame(i)
    n_frame += 1
    for sum_atom, atom in zip(sum_soup.atoms(), trj.soup.atoms()):
      sum_atom.pos += atom.pos
  
  for a in sum_soup.atoms():
    x = a.pos.x / float(n_frame)
    y = a.pos.y / float(n_frame)
    z = a.pos.z / float(n_frame)
    v3.set_vector(a.pos, x, y, z)
  
  trj.soup.write_pdb(out_name)


def get_n_frame_per_ps(name):
  config = name + ".config"
  if not os.path.isfile(config):
    config = "pulse0/" + name + ".config"
    if not os.path.isfile(config):
      raise ValueError("Can't find .config file")
  parms = util.read_dict(config)
  if 'n_step_per_snapshot' in parms:
    n_step_per_snapshot = parms['n_step_per_snapshot']
  elif 'steps_per_snapshot' in parms:
    n_step_per_snapshot = parms['steps_per_snapshot']
  return 1000 / n_step_per_snapshot
  
  
def make_pdb_from_traj(md_name, ps, out_pdb, res_vals=None):
  n_frame_per_ps = get_n_frame_per_ps(md_name)
  traj = open_trajectory(md_name)
  if ps == 0:
    i_frame = 0
  else:
    i_frame = int(ps * n_frame_per_ps) - 1
  traj.load_frame(i_frame)
  
  if res_vals is None:
    traj.soup.write_pdb(out_pdb)
  else:
    traj.soup.load_residue_bfactors(res_vals)
    traj.soup.write_pdb(out_pdb)
  del traj

  
def get_non_solvent_residues(soup):
  residues = []
  for r in soup.residues():
    if r.type not in data.solvent_res_types:
      residues.append(r)
  return residues
  
  
class ProcessVariable:
  def set(self, trj, n_frame_per_ps, fname, ref_pdb):
    self.fname = trj.name + '.' + fname
    self.trj = trj
    self.n_frame_per_ps = n_frame_per_ps
    self.n_ps = self.trj.n_frame / self.n_frame_per_ps
    self.residues = get_non_solvent_residues(self.trj.soup)
    self.n_residue = len(self.residues)
    self.cumul_res_averages = [0.0 for i in range(self.n_residue)]
    self.res_averages_per_ps = []
    self.atoms = []
    for r in self.residues:
      self.atoms.extend(r.atoms())
    self.n_atom = len(self.atoms)

  def is_done(self):
    return os.path.isfile(self.fname)
  
  def open(self):
    self.file = open(self.fname, 'w')

  def process_frame_on_ps(self):
    if self.trj.i_frame == 0:
      n_frame = 1
    else:
      n_frame = self.n_frame_per_ps
    # save averaged sc values in array for later
    for i in range(self.n_residue):
      self.cumul_res_averages[i] /= float(n_frame)
    self.res_averages_per_ps.append(
        copy.deepcopy(self.cumul_res_averages))
    
    # clear for next ps
    if self.trj.i_frame > 0:
      for j in range(self.n_residue):
        self.cumul_res_averages[j] = 0.0

  def extract_vals_to_atom(self):
    # write to self.trj.res(i).atom.val
    pass
  
  def process_frame(self):
    time = self.trj.i_frame / float(self.n_frame_per_ps)
    self.extract_vals_to_atom()
    res_averages = [0.0 for i in range(self.n_residue)]
    for i in range(self.n_residue):
      atoms = self.residues[i].atoms()
      res_sum = sum([atom.val for atom in atoms])
      res_averages[i] = res_sum / float(len(atoms))
    for i in range(self.n_residue):
      self.cumul_res_averages[i] += res_averages[i]
    # write res_averages to file
    res_str = " ".join([str(val) for val in res_averages])
    self.file.write("%f %s\n" % (time, res_str))
    
  def print_res_averages_per_ps(self):
    f = open('%s_ave' % self.fname, "w")
    for i in range(self.n_residue):
      vals_by_time = [res_vals[i] 
                      for res_vals in self.res_averages_per_ps]
      str_list = [str(v) for v in vals_by_time]
      f.write(' '.join(str_list) + "\n")
    f.close()

  def close(self):
    self.file.close()
    self.print_res_averages_per_ps()


class ProcessRmsd(ProcessVariable):
  def set(self, trj, n_frame_per_ps, ref_pdb):
    ProcessVariable.set(self, trj, n_frame_per_ps, "rmsd", ref_pdb)
    if ref_pdb:
      self.ref_soup = pdbatoms.Polymer(ref_pdb)
    else:
      self.ref_soup = self.trj.soup.copy()
    self.ref_soup.write_pdb("ref.pdb")
    self.ref_residues = get_non_solvent_residues(self.ref_soup)
    self.ref_atoms = []
    for r in self.ref_residues:
      self.ref_atoms.extend(r.atoms())
  
  def extract_vals_to_atom(self):
    for atom, ref_atom in zip(self.atoms, self.ref_atoms):
      atom.val = v3.distance(atom.pos, ref_atom.pos)
    
  
class ProcessKineticEnergy(ProcessVariable):
  def set(self, trj, n_frame_per_ps, ref_pdb):
    ProcessVariable.set(self, trj, n_frame_per_ps, "kin", ref_pdb)
  
  def extract_vals_to_atom(self):
    for atom in self.atoms:
      vel = v3.mag(atom.vel)
      atom.val = 0.5 * atom.mass * vel * vel


class ProcessCaRmsd(ProcessVariable):
  def set(self, trj, n_frame_per_ps, ref_pdb):
    ProcessVariable.set(self, trj, n_frame_per_ps,
                        "ca_rmsd", ref_pdb)
    if ref_pdb:
      self.ref_soup = pdbatoms.Polymer(ref_pdb)
    else:
      self.ref_soup = self.trj.soup.copy()
    self.ref_residues = get_non_solvent_residues(self.ref_soup)
    self.ref_atoms = []
    for r in self.ref_residues:
      self.ref_atoms.extend(r.atoms())

  def process_frame(self):
    time = self.trj.i_frame / float(self.n_frame_per_ps)

    res_averages = [0.0 for i in range(self.n_residue)]
    for i in range(self.n_residue):
      if self.residues[i].has_atom('CA'):
        res_averages[i] = v3.distance(
            self.residues[i].atom('CA').pos,
            self.ref_residues[i].atom('CA').pos)
    res_str = " ".join([str(val) for val in res_averages])
    self.file.write("%f %s\n" % (time, res_str))

    for i in range(self.n_residue):
      self.cumul_res_averages[i] += res_averages[i]


class ProcessTotalKineticEnergy(ProcessVariable):
  def set(self, trj, n_frame_per_ps, ref_pdb):
    ProcessVariable.set(self, trj, n_frame_per_ps, 
                        "total_kin", ref_pdb)
    self.n_atom = len(trj.soup.atoms())

  def process_frame(self):
    time = self.trj.i_frame / float(self.n_frame_per_ps)
    energy = 0.0
    for atom in self.atoms:
      vel = v3.mag(atom.vel)
      energy += 0.5 * atom.mass * vel * vel
    self.file.write("%f %f %f\n" % \
                    (time, energy, energy / self.n_atom))

  def print_res_averages_per_ps(self):
    pass


def process_trajectory(
    md, n_frame_per_ps, in_processes, ref_pdb=None):
  """"
  Looks for (name.psf, name.dcd, name.vel.dcd) or
  (name.top, name.trj, name.vel) to process trajectories
  """
  trj =  open_trajectory(md)
  processes = in_processes  
  for process in processes:
    process.set(trj, n_frame_per_ps, ref_pdb)
  processes = [process for process in processes \
               if not process.is_done()]
  if processes:
    for process in processes:
      process.open()
    for i in range(trj.n_frame):
      trj.load_frame(i)
      for process in processes:
        process.process_frame()
      if (i+1) % n_frame_per_ps == 0 or i == 0:
        for process in processes:
          process.process_frame_on_ps()
    for process in processes:
      process.close()


def process_md(md, ref_pdb=None):
  config = md + '.config'
  if not os.path.isfile(config):
    n_frame_per_ps = 50 # default if not given
  else:
    n_frame_per_ps = get_n_frame_per_ps(md)
  processes = [ProcessKineticEnergy(), 
               ProcessRmsd(), 
               ProcessCaRmsd(),
               ProcessTotalKineticEnergy()]
  process_trajectory(md, n_frame_per_ps, processes, ref_pdb)

