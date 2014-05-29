

__doc__ = """
This is an example of a system equilibration to 10 K
and then applying a RIP rotational perturbation on a 
residue using the pdbremix library.
"""

import os

from pdbremix import pdbatoms
from pdbremix import simulate
from pdbremix import force
from pdbremix import util
from pdbremix import pdbtext
from pdbremix import fetch
from pdbremix import data


def make_protein_restraint_pdb(in_basename, pdb):
  soup = simulate.soup_from_restart_files(in_basename)
  for a in soup.atoms():
    a.bfactor = 0.0
    if a.res_type not in data.solvent_res_types:
      a.bfactor = 1.0
  soup.write_pdb(pdb)


sim_dir = "equilibrate"
pdb_code = '1cph'
ff = 'GROMACS4.5'
pdb_code = '2evq'
ff = 'AMBER11-GBSA'
ff = 'GROMACS4.5'
i_residue = 2

util.goto_dir(sim_dir)

# Get the PDB files from the website
fetch.get_pdbs_with_http(pdb_code)

# Format the PDB file to give one single conformation 
clean_pdb = '2evq.clean.pdb'
pdbtext.clean_pdb(pdb_code+'.pdb', clean_pdb)

# Generate restart files from PDB
simulate.pdb_to_top_and_crds(ff, clean_pdb, 'sim')

# Make a protein positional restrain_protein.pdb file
make_protein_restraint_pdb('sim', 'restrain_protein.pdb')

# minimize system (mostly water) with protein positions fixed
simulate.minimize(
    ff, 'sim', 'minwater', restraint_pdb='restrain_protein.pdb')

# heat only water to 300K, hold protein fixed
simulate.langevin_thermometer(
    ff, 'minwater', 5000, 300, 'heatwater',  50, 
    restraint_pdb='restrain_protein.pdb')

# cool water back down 10K, hold protein fixed
simulate.langevin_thermometer(
    ff, 'heatwater', 5000, 10, 'coolwater',  50, 
    restraint_pdb='restrain_protein.pdb')

# equilibrate entire system to 10K 
simulate.langevin_thermometer(
    ff, 'coolwater', 5000, 10, 'heat',  50)

# let the system relax without thermometer
simulate.constant_energy(
    ff, 'heat', 5000, 'const',  50)

# then reequilibrate to 10K
simulate.langevin_thermometer(
    ff, 'const', 5000, 10, 'reheat',  50)

# we are ready to apply RIP
pulse_fn = force.make_rip_fn(i_residue, 100)
simulate.pulse(
    ff, 'reheat', 'rip', 5000, pulse_fn, 100)

# combine all sims into one long trajectory for viewing
simulate.merge_trajectories(
    ff, 'equil', ['heatwater', 'coolwater', 'heat', 
                  'const', 'reheat', 'rip'])
