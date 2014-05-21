import os

from pdbremix import pdbatoms
from pdbremix import simulate
from pdbremix import force
from pdbremix import util
from pdbremix import pdbtext
from pdbremix import trajectory
from pdbremix import fetch
from pdbremix import v3


def test_prepare_for_md(params):
  util.goto_dir(sim_dir)

  print os.getcwd()
  pdb = params['pdb']

  fetch.get_pdbs_with_http(pdb)
  full_pdb = pdb + '.pdb'
  clean_pdb = os.path.abspath(pdb) + '.clean.pdb'
  pdbtext.clean_pdb(full_pdb, clean_pdb)

  print "> Generating topologies"
  top, crds = simulate.pdb_to_top_and_crds(
      params['ff'], clean_pdb, 'sim')
  util.goto_dir('min')

  print "> Minimizing structure"
  top, crds, vels = simulate.get_restart_files('../sim')
  simulate.minimize(params['ff'], top, crds, 'min')


def test_basic_md_merge(params):
  util.goto_dir(params['sim_dir'])

  print "> Fixed-temp1"
  util.goto_dir('md1')
  top, crds, vels = simulate.get_restart_files('../min/min')
  simulate.langevin_thermometer(
      params['ff'], top, crds, vels, 1000, 300, 'md', 10)

  print "> Fixed-temp2"
  util.goto_dir('../md2')
  top, crds, vels = simulate.get_restart_files('../md1/md')
  simulate.langevin_thermometer(
      params['ff'], top, crds, vels, 1000, 300, 'md', 10)

  print "> Merge temp1-temp2"
  util.goto_dir('../md_merge')
  simulate.merge_simulations(
      params['ff'], 'md', ['../md1', '../md2'])


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
    # make it bounce
    if self.last_rot_vel*rot_vel < 0.0:
      self.rot_vel = -self.rot_vel
    force.add_rotational_velocity(
        move_atoms, self.rot_vel - rot_vel, axis, anchor1)


def test_user_defined_pulse(params):
  util.goto_dir(params['sim_dir'])
  util.goto_dir('user_pulse')
  print "> Pulse with user-defined pulse_fn"
  rotate_loop = RotateLoop(39, 44, -0.7)
  pulse_fn = lambda soup: rotate_loop.apply(soup)
  print "> Pulse with terminii push"
  simulate.pulse(
      params['ff'], '../md_merge/md', 'md', 5000, pulse_fn, 100)


def test_rip(params):
  util.goto_dir(params['sim_dir'])
  print "> Pulse with RIP residue rotation"
  util.goto_dir('rip')
  pulse_fn = force.make_rip_fn(params['i_residue'], 300)
  simulate.pulse(
      params['ff'], '../md_merge/md', 'md', 2000, pulse_fn, 100)


def test_puff(params):
  util.goto_dir(params['sim_dir'])
  util.goto_dir('puff')

  # make push function
  md = '../md_merge/md'
  n = len(simulate.soup_from_restart_files(md).residues())
  pulse_fn = force.make_puff_fn(
      [0, 1, 2], [n-3, n-2, n-1], 10.0, 0.1, 300)

  print "> Pulse with terminii push"
  simulate.pulse(params['ff'], md, 'md', 2000, pulse_fn, 100)


def test_restraint(params):
  util.goto_dir(params['sim_dir'])

  util.goto_dir('restraint')

  md = '../md_merge/md'

  # make restraint_pdb
  restraint_pdb = 'restraint.pdb'
  soup = simulate.soup_from_restart_files(md)
  for atom in soup.residue(params['i_residue']).atoms():
    atom.bfactor = 1.0
  soup.write_pdb(restraint_pdb)

  print "> Equil with restraints on "
  top, crds, vels = simulate.get_restart_files(md)
  simulate.langevin_thermometer(
      params['ff'], top, crds, vels, 1000, 300, 'md',  10, 
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


def test_all(params):
  test_prepare_for_md(params)
  test_basic_md_merge(params)
  test_user_defined_pulse(params)
  test_rip(params)
  test_puff(params)
  test_restraint(params)
  test_traj_analysis(params)


if __name__ == "__main__":

  # ff = 'AMBER11-GBSA'
  # ff = 'AMBER11'
  # ff = 'GROMACS4.5'
  # ff = 'NAMD2.8'

  # # insulin
  # pdb = '1cph' 
  # i_residue = 18

  # # beta-hairpin turn
  # pdb = '2evq'
  # i_residue = 2 

  params = {
    'ff': 'AMBER11-GBSA',
    'pdb': '1cph',
    'i_residue': 2,
    'save_dir': os.getcwd()
  }
  sim_dir = '%s/%s' % (params['ff'], params['pdb'])
  params['sim_dir'] = os.path.abspath(sim_dir)

  # util.clean_fname(sim_dir)

  test_all(params)

