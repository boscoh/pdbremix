

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
from pdbremix import trajectory
from pdbremix import fetch
from pdbremix import data

sim_dir = "equilibrate"
pdb = '2evq'
pdb = '1cph'
ff = 'GROMACS4.5'

util.goto_dir(sim_dir)


# Get the PDB files from the website
fetch.get_pdbs_with_http(pdb)

# Format the PDB file to give one single conformation 
clean_pdb = '2evq.clean.pdb'
pdbtext.clean_pdb(pdb+'.pdb', clean_pdb)

# Generate restart files from PDB
top, crds = simulate.pdb_to_top_and_crds(ff, clean_pdb, 'sim')

# Make a protein positional restraint.pdb file
soup = simulate.soup_from_restart_files('sim')
for a in soup.atoms():
  a.bfactor = 0.0
  if a.res_type not in data.solvent_res_types:
    a.bfactor = 1.0
soup.write_pdb('restraint.pdb')

# minimize system (mostly water) with protein positions fixed
util.goto_dir('min')
top, crds, vels = simulate.get_restart_files('../sim')
simulate.minimize(
    ff, top, crds, 'min', 
    restraint_pdb='../restraint.pdb')

# heat only water to 300K, hold protein fixed
util.goto_dir('../heatwater')
top, crds, vels = simulate.get_restart_files('../min/min')
simulate.langevin_thermometer(
    ff, top, crds, vels, 5000, 300, 'md',  50, 
    restraint_pdb='../restraint.pdb')

# cool water back down 10K, hold protein fixed
util.goto_dir('../coolwater')
top, crds, vels = simulate.get_restart_files('../heatwater/md')
simulate.langevin_thermometer(
    ff, top, crds, vels, 5000, 10, 'md',  50, 
    restraint_pdb='../restraint.pdb')

# equilibrate entire system to 10K 
util.goto_dir('../heat')
top, crds, vels = simulate.get_restart_files('../coolwater/md')
simulate.langevin_thermometer(
    ff, top, crds, vels, 5000, 10, 'md',  50)

# let the system relax without thermometer
util.goto_dir('../const')
top, crds, vels = simulate.get_restart_files('../heat/md')
simulate.constant_energy(
    ff, top, crds, vels, 5000, 'md',  50)

# then reequilibrate to 10K
util.goto_dir('../reheat')
top, crds, vels = simulate.get_restart_files('../const/md')
simulate.langevin_thermometer(
    ff, top, crds, vels, 5000, 10, 'md',  50)

# we are ready to apply RIP
util.goto_dir('../rip')
simulate.pulse(
    ff, '../reheat/md', 'md', 5000, 
    force.make_rip_fn(18, 100), 100)

# splice all the simulations together into one trajectory
util.goto_dir('../equil')
simulate.merge_simulations(
    ff, 'md', 
    ['../heatwater', '../coolwater', '../heat', 
     '../const', '../reheat', '../rip'])
