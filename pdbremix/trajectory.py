# encoding: utf-8

__doc__ = """
This module provides a way to abstract the reading of 
Molecular-Dynamics (MD) trajectories. 

Simulation trajectories are read into a Trajectory object
that provides a Soup object, which is updated everytime a
frame is requested. All operations on Soup are available,
including the writing of PDB files, and extracting fragments
for further interrogation.

Trajectory object:
  Attributes:
    basename
    n_frame
    i_frame
    soup
  Methods:
    load_frame

Currently supported trajectories:
- AMBER: .top .trj .vel.trj
- GROMACS: .top .gro .trr
- NAMD: .psf .gro .dcd .vel.dcd
"""

import os
import copy

import util
import v3
import pdbatoms
import data

import namd
import amber
import gromacs



def open_trajectory(basename):
  """
  Returns a Trajectory object by guessing the type from basename.
  """
  if os.path.isfile(basename + '.trr'):
    print "Detected GROMACS trajectory"
    return gromacs.Trajectory(basename)
  if os.path.isfile(basename + '.psf'):
    print "Detected NAMD trajectory"
    return namd.Trajectory(basename)
  if os.path.isfile(basename + '.top'):
    print "Detected AMBER trajectory"
    return amber.Trajectory(basename)
  raise ValueError("No trajectories found at all!")


def make_equilibrium_pdb(trajectory, pdb, start=None, end=None):
  """
  Calculates the equilibrium structure between start and end
  frames. Default: last half of trajectory.
  """
  if start is None:
    start = trajectory.n_frame // 2
  if end is None:
    end = trajectory.n_frame
  
  if start < 0 or start >= trajectory.n_frame:
    raise IndexError("Start frame out of range")
  if end < 0 or end > trajectory.n_frame+1:
    raise IndexError("End frame out of range")
  
  # initialize sum object
  sum_soup = trajectory.soup.copy()
  for a in sum_soup.atoms():
    v3.set_vector(a.pos, 0, 0, 0)
  
  n_frame = 0
  for i in range(start, end):
    trajectory.load_frame(i)
    n_frame += 1
    for sum_atom, atom in zip(
        sum_soup.atoms(), trajectory.soup.atoms()):
      sum_atom.pos += atom.pos
  
  for a in sum_soup.atoms():
    x = a.pos.x / float(n_frame)
    y = a.pos.y / float(n_frame)
    z = a.pos.z / float(n_frame)
    v3.set_vector(a.pos, x, y, z)
  
  trajectory.soup.write_pdb(pdb)


def make_pdb_from_trajectory(
    basename, i_frame, out_pdb, res_bfactors=None):
  """
  Make a PDB from the ith frame of a trajectory described by
  basename. Optional residue bfactors are loaded into residues.
  """
  trajectory = open_trajectory(basename)
  trajectory.load_frame(i_frame)
  if res_vals is not None:
    trajectory.soup.load_residue_bfactors(res_bfactors)
  trajectory.soup.write_pdb(out_pdb)
  del trajectory



def get_non_solvent_residues(soup):  
  residues = []
  for r in soup.residues():
    if r.type not in data.solvent_res_types:
      residues.append(r)
  return residues



class TrajectoryAnalyzer(object):
  """
  Abstract strategy object to analyze the frames of a trajectory
  and write the results to a text file.
  """

  def __init__(self, trj, n_frame_per_ps, var_name, ref_pdb):
    self.trj = trj
    self.n_frame_per_ps = n_frame_per_ps

    # collect non-solvent residues
    self.residues = get_non_solvent_residues(self.trj.soup)
    self.n_residue = len(self.residues)
    self.cumul_res_averages = [0.0 for i in range(self.n_residue)]
    self.res_averages_per_ps = []

    self.atoms = self.trj.soup.atoms()
    self.n_atom = len(self.atoms)

    self.fname = trj.basename + '.' + var_name + '.per_frame'
    self.file = open(self.fname, 'w')

  def process_frame_on_ps(self):
    if self.trj.i_frame == 0:
      n_frame = 1
    else:
      n_frame = self.n_frame_per_ps

    # save averaged res values in array for later
    for i in range(self.n_residue):
      self.cumul_res_averages[i] /= float(n_frame)
    self.res_averages_per_ps.append(
        copy.deepcopy(self.cumul_res_averages))
    
    # clear for next ps
    if self.trj.i_frame > 0:
      for j in range(self.n_residue):
        self.cumul_res_averages[j] = 0.0

  def extract_vals_to_atom(self):
    """
    To be overriden by desccendent classes. Results
    are to be written to atom.val for each atom.
    """
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
    f = open(self.fname.replace('frame', 'ps'), "w")
    for i in range(self.n_residue):
      vals_by_time = [res_vals[i] 
                      for res_vals in self.res_averages_per_ps]
      str_list = [str(v) for v in vals_by_time]
      f.write(' '.join(str_list) + "\n")
    f.close()

  def close(self):
    self.file.close()
    self.print_res_averages_per_ps()



class RmsdAnalyzer(TrajectoryAnalyzer):
  """
  TrajectoryAnalyzer object to calculate the RMSD of residues.
  """
  def __init__(self, trj, n_frame_per_ps, ref_pdb):
    TrajectoryAnalyzer.__init__(self, trj, n_frame_per_ps, "rmsd", ref_pdb)
    if ref_pdb:
      self.ref_soup = pdbatoms.Soup(ref_pdb)
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
    

  
class CaRmsdAnalyzer(TrajectoryAnalyzer):
  """
  TrajectoryAnalyzer object to calculate the C-alpha RMSD.
  """
  def __init__(self, trj, n_frame_per_ps, ref_pdb):
    TrajectoryAnalyzer.__init__(
        self, trj, n_frame_per_ps, "carmsd", ref_pdb)
    if ref_pdb:
      self.ref_soup = pdbatoms.Soup(ref_pdb)
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



class KineticEnergyAnalyzer(TrajectoryAnalyzer):
  """
  TrajectoryAnalyzer object to calculate kinetic energy 
  of residues.
  """
  def __init__(self, trj, n_frame_per_ps, ref_pdb):
    TrajectoryAnalyzer.__init__(
      self, trj, n_frame_per_ps, "kin", ref_pdb)
  
  def extract_vals_to_atom(self):
    for atom in self.atoms:
      vel = v3.mag(atom.vel)
      atom.val = 0.5 * atom.mass * vel * vel



class TotalKineticEnergyAnalyzer(TrajectoryAnalyzer):
  """
  Strategy object to calculate the total kinetic energy.
  """
  def __init__(self, trj, n_frame_per_ps, ref_pdb):
    TrajectoryAnalyzer.__init__(
        self, trj, n_frame_per_ps, "total_kin", ref_pdb)

  def process_frame(self):
    time = self.trj.i_frame / float(self.n_frame_per_ps)
    energy = 0.0
    for atom in self.atoms:
      vel = v3.mag(atom.vel)
      energy += 0.5 * atom.mass * vel * vel
    self.file.write(
        "%f %f %f\n" % (time, energy, energy / self.n_atom))

  def print_res_averages_per_ps(self):
    """
    Override, not needed. 
    """
    pass



def guess_n_frame_per_ps(basename):
  """
  Returns the n_frame_per_ps of a trajectory by reading any
  .config files that would have been generated using simualte.py.
  """
  config = basename + ".config"
  try:
    params = util.read_dict(config)
    # assuming 1fs time step
    n_step_per_ps = 1000 
    if 'n_step_per_snapshot' in params:
      n_step_per_snapshot = params['n_step_per_snapshot']
    n_frame_per_ps = n_step_per_ps / n_step_per_snapshot
  except:
    n_frame_per_ps = 50
  return n_frame_per_ps
  
  

def analyze_trajectory(
    basename, n_frame_per_ps=None, 
    analyzer_classes=None, ref_pdb=None):
  """"
  Calculates trajectory variables for the residues in the 
  system using the TrajectoryAnalyzer straegy objects in 
  analyzers.
  """
  if n_frame_per_ps is None:
    n_frame_per_ps = guess_n_frame_per_ps(basename)

  if analyzer_classes is None:
    analyzer_classes = [
        KineticEnergyAnalyzer, 
        RmsdAnalyzer, 
        CaRmsdAnalyzer,
        TotalKineticEnergyAnalyzer]

  trj = open_trajectory(basename)

  # Instantiate analyzers
  analyzers = [a(trj, n_frame_per_ps, ref_pdb) \
               for a in analyzer_classes]

  # Go through frame by frame
  for i in range(trj.n_frame):
    trj.load_frame(i)
    for analyzer in analyzers:
      analyzer.process_frame()
    if (i+1) % n_frame_per_ps == 0 or i == 0:
      for analyzer in analyzers:
        analyzer.process_frame_on_ps()

  # Close analyzers
  for analyzer in analyzers:
    analyzer.close()



