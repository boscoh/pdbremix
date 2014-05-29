
__doc__ = """

Example scripts to run MD simulations in PDBREMIX.

It runs through a bunch of different simulations options, with
user-defined pulsing functions,  and user-defined trajectory
analysis. Script will fail on exception if there are any
problems.

"""

# Setup parameters - force_field and structure

# Alternate SMALL structures
# insulin: 2-chains & 3 disuflide bonds
#    'pdb': '1cph', 
#    'i_residue': 18,
#    'i_loop_start': 39, 
#    'i_loop_end': 44
# beta-hairpin: 
#    'pdb': '2evq', 
#    'i_residue': 2, 

# Available force-fields:
#   AMBER11-GBSA, AMBER11, GROMACS4.5, NAMD2.8

params = {
  'ff': 'AMBER11-GBSA',
  'pdb': '2evq',
  'i_residue': 2,
  'i_loop_start': 3, 
  'i_loop_end': 7
}


# 'ff': 'GROMACS4.5',

params = {
  'ff': 'AMBER11-GBSA',
  'pdb': '1cph',
  'i_residue': 18,
  'i_loop_start': 39, 
  'i_loop_end': 44
}

params = {
  'ff': 'AMBER11-GBSA',
  'pdb': '2evq',
  'i_residue': 2,
  'i_loop_start': 3, 
  'i_loop_end': 7
}

import os

from pdbremix import pdbatoms
from pdbremix import simulate
from pdbremix import force
from pdbremix import util
from pdbremix import pdbtext
from pdbremix import trajectory
from pdbremix import fetch
from pdbremix import data
from pdbremix import v3


def make_restraint_pdb(in_md, residue_indices, out_pdb, is_backbone_only=True):
  soup = simulate.soup_from_restart_files(in_md)
  for i_res in residue_indices:
    for a in soup.residue(i_res).atoms():
      if is_backbone_only and a.type not in data.backbone_atoms:
        continue
      a.bfactor = 1.0
  soup.write_pdb(out_pdb)


def get_n_protein_residues_from_restart(in_md):
  n = 0
  soup = simulate.soup_from_restart_files(in_md)
  for res in soup.residues():
    if res not in data.solvent_res_types:
      n += 1
  return n


def test_prepare_for_md(params):
  util.goto_dir(sim_dir)

  print os.getcwd()
  pdb = params['pdb']

  fetch.get_pdbs_with_http(pdb)
  full_pdb = pdb + '.pdb'
  clean_pdb = os.path.abspath(pdb) + '.clean.pdb'
  pdbtext.clean_pdb(full_pdb, clean_pdb)

  print "> Generating topologies"
  simulate.pdb_to_top_and_crds(params['ff'], clean_pdb, 'sim')

  print "> Minimizing structure"
  util.goto_dir('min')
  simulate.minimize(params['ff'], '../sim', 'min')


def test_basic_md_merge(params):
  util.goto_dir(params['sim_dir'])

  print "> Fixed-temp1"
  util.goto_dir('md1')
  simulate.langevin_thermometer(
      params['ff'], '../min/min', 1000, 300, 'md', 10)

  print "> Fixed-temp2"
  util.goto_dir('../md2')
  simulate.langevin_thermometer(
      params['ff'], '../md1/md', 1000, 300, 'md', 10)

  print "> Merge temp1-temp2"
  util.goto_dir('../md_merge')
  simulate.merge_trajectories(
      params['ff'], 'md', ['../md1/md', '../md2/md'])


class RotateLoop:
  def __init__(self, i, j, rot_vel):
    self.i, self.j = i, j
    self.last_rot_vel = 0.0
    self.rot_vel = rot_vel
  def apply(self, soup):
    anchor1 = soup.residue(self.i).atom('CA').pos
    anchor2 = soup.residue(self.j).atom('CA').pos
    axis = anchor1 - anchor2
    move_atoms = []
    for k in range(self.i+1, self.j):
      move_atoms.extend(soup.residue(k).atoms())
    rot_vel = force.weighted_rotational_velocity(
        move_atoms, axis, anchor1)
    force.add_rotational_velocity(
        move_atoms, self.rot_vel - rot_vel, axis, anchor1)


def test_user_defined_pulse(params):
  util.goto_dir(params['sim_dir'])
  util.goto_dir('user_pulse')
  i, j = params['i_loop_start'], params['i_loop_end']
  rotate_loop = RotateLoop(i, j, 0.7)
  print "> make restraint_pdb"
  restraint_pdb = os.path.abspath('md.restraint.pdb')
  in_md = '../md_merge/md'
  n = get_n_protein_residues_from_restart(in_md)
  residue_indices = range(0, i) + range(j, n)
  make_restraint_pdb(in_md, residue_indices, restraint_pdb)
  pulse_fn = lambda soup: rotate_loop.apply(soup)
  print "> Pulse with user-defined pulse_fn"
  simulate.pulse(
      params['ff'], '../md_merge/md', 'md', 5000, pulse_fn, 100,
      restraint_pdb=restraint_pdb)


def test_rip(params):
  util.goto_dir(params['sim_dir'])
  print "> Pulse with RIP residue rotation"
  util.goto_dir('rip')
  simulate.pulse(
      params['ff'], '../md_merge/md', 'md', 2000, 
      force.make_rip_fn(params['i_residue'], 300), 100)


def test_puff(params):
  util.goto_dir(params['sim_dir'])
  util.goto_dir('puff')
  print "> Pulse with terminii push"
  md = '../md_merge/md'
  n = get_n_protein_residues_from_restart(md)
  # make push function
  pulse_fn = force.make_puff_fn(
      [0, 1, 2], [n-3, n-2, n-1], 10.0, 0.1, 300)
  simulate.pulse(params['ff'], md, 'md', 2000, pulse_fn, 100)


def test_restraint(params):
  print "> Equil with restraints on "
  util.goto_dir(params['sim_dir'])
  util.goto_dir('restraint')
  in_md = '../md_merge/md'
  restraint_pdb = os.path.abspath('md.restraint.pdb')
  make_restraint_pdb(
      in_md, [params['i_residue']], restraint_pdb, is_backbone_only=False)
  simulate.langevin_thermometer(
      params['ff'], in_md, 1000, 300, 'md',  10, 
      restraint_pdb=restraint_pdb)


class VelocityAnalyzer(trajectory.TrajectoryAnalyzer):
  var_name = 'velocity'
  def calculate_results(self):
    atoms = self.soup.atoms()
    n = len(atoms)
    vel = v3.vector()
    for a in self.soup.atoms():
      vel += a.vel
    vel = v3.scale(vel, float(n))
    return [float("%02.f" % v) for v in vel]


def test_traj_analysis(params):
  util.goto_dir(params['sim_dir'])
  util.goto_dir('md1')
  trajectory.analyze_trajectory('md')
  util.check_output('md.rmsd.per_frame')
  util.check_output('md.rmsd.per_ps')

  trajectory.analyze_trajectory(
      'md', analyzer_classes=[VelocityAnalyzer])
  util.check_output('md.velocity.per_frame')
  util.check_output('md.velocity.per_ps')



sim_dir = '%s/%s' % (params['ff'], params['pdb'])
params.update({
  'sim_dir': os.path.abspath(sim_dir),
  'save_dir': os.getcwd(),
})
test_prepare_for_md(params)
test_basic_md_merge(params)
test_user_defined_pulse(params)
test_rip(params)
test_puff(params)
test_restraint(params)
test_traj_analysis(params)
