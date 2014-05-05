import os

from pdbremix import pdbatoms
from pdbremix import simulate
from pdbremix import force
from pdbremix import util


ff = 'GROMACS4.5'
ff = 'NAMD2.8'
ff = 'AMBER11-GBSA'
ff = 'AMBER11'
ff = 'GROMACS4.5'
pdb = 'pdb/1cph.pdb'
pdb = 'pdb/hairpin.pdb'

pdb = os.path.abspath(pdb)
name = os.path.splitext(os.path.basename(pdb))[0]
sim_dir = 'md/%s/%s' % (ff, name)

# util.clean_fname(sim_dir)

save_dir = os.getcwd()


def prepare_for_md():
  util.goto_dir(sim_dir)
  top, crds = simulate.pdb_to_top_and_crds(ff, pdb, 'sim')
  util.goto_dir('min')
  top, crds, vels = simulate.get_restart_files('../sim')
  simulate.minimize(ff, top, crds, 'min')
  util.goto_dir(save_dir)


def test_prepare_for_md():
  prepare_for_md()


def test_basic_md_merge():
  prepare_for_md()
  util.goto_dir(sim_dir)
  util.goto_dir('md1')
  top, crds, vels = simulate.get_restart_files('../min/min')
  simulate.langevin(ff, top, crds, vels, 1000, 300, 'md', 10)
  util.goto_dir('..')
  util.goto_dir('md2')
  top, crds, vels = simulate.get_restart_files('../md1/md')
  simulate.langevin(ff, top, crds, vels, 1000, 300, 'md', 10)
  util.goto_dir('..')
  util.goto_dir('md_merge')
  simulate.merge_simulations(ff, 'md', ['../md1', '../md2'])
  util.goto_dir('..')
  util.goto_dir(save_dir)


def test_rip():
  prepare_for_md()
  util.goto_dir(sim_dir)
  util.goto_dir('rip')
  md = '../md_merge/md'
  top, crds, vels = simulate.get_restart_files(md)
  soup = simulate.soup_from_restart_files(top, crds, vels)
  i_residue = 2
  res = soup.residue(i_residue)
  pulse_fn = force.make_rip_fn(i_residue, 300)
  simulate.pulse(ff, md, 'md', 2000, pulse_fn, 100)
  util.goto_dir(save_dir)


def test_puff():
  prepare_for_md()
  util.goto_dir(sim_dir)
  util.goto_dir('puff')
  md = '../md_merge/md'
  top, crds, vels = simulate.get_restart_files(md)
  n = len(simulate.soup_from_restart_files(top, crds, vels).residues())
  pulse_fn = force.make_puff_fn([0, 1, 2], [n-3, n-2, n-1], 10.0, 0.1, 300)
  simulate.pulse(ff, md, 'md', 2000, pulse_fn, 100)
  util.goto_dir(save_dir)

def test_restraint():
  prepare_for_md()
  util.goto_dir(sim_dir)
  util.goto_dir('restraint')
  top, crds, vels = simulate.get_restart_files('../min/min')
  restraint_pdb = '../sim.restraint.pdb'
  simulate.convert_restart_to_pdb('../min/min', restraint_pdb)
  soup = pdbatoms.Polymer(restraint_pdb)
  for atom in soup.residue(2).atoms():
    atom.bfactor = 1.0
  soup.write_pdb(restraint_pdb)
  simulate.langevin(
      ff, top, crds, vels, 1000, 300, 'md',  10, restraint_pdb=restraint_pdb)
  util.goto_dir(save_dir)


if __name__ == "__main__":
  test_basic_md_merge()
  test_rip()
  test_puff()
  test_restraint()


