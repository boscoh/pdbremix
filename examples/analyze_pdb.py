
__doc__ = """

Example scripts to analyze PDB in PDBREMIX.

This goes through a number of analysis that are reasonably
useful. It demonstrates how the PDBREMIX data-structures 
interfaces with Python, using idiomatic Python approaches
such as lambda, list-comprehensions, re module, and
a generally more functional-approach.

"""

import os
import re

from pdbremix import v3
from pdbremix import pdbatoms
from pdbremix import simulate
from pdbremix import force
from pdbremix import util
from pdbremix import pdbtext
from pdbremix import trajectory
from pdbremix import fetch
from pdbremix import asa
from pdbremix import data
from pdbremix import rmsd


util.goto_dir('pdb')


# load the structures from the website
fetch.get_pdbs_with_http('1be9', '1rzx')
soup = pdbatoms.Soup('1be9.pdb')
soup2 = pdbatoms.Soup('1rzx.pdb')



# the structures loaded 1be9 and 1rzx are both PDZ
# domains. However, they don't quite map to each other
# This shows how you might select homologous residues
# and then extract the coordinates for RMSD calculation

def split_pairs(indices):
  for i in range(0, len(indices), 2):
    yield indices[i:i+2]

def get_crds(res_tags, soup):
  crds = []
  i_residues = [soup.get_i_residue(r) for r in res_tags]
  for j, k in split_pairs(i_residues):
    for i in range(j, k+1):
      crds.append(soup.residue(i).atom('CA').pos)
  center = v3.get_center(crds)
  crds = [c-center for c in crds]    
  return crds

res_tags1 = "A:311 A:318 A:322 A:329 A:335 A:391".split()
res_tags2 = "A:158 A:165 A:171 A:178 A:194 A:250".split()
rmsd, rot = rmsd.pyqcprot_rmsd_rot(
    get_crds(res_tags1, soup), get_crds(res_tags2, soup2))
print "RMSD:"
print rmsd



# select polar/charged residues using regular expresions
# and list comprehensions
residues = filter(
    lambda r: re.match('(LY.|AR.|GL.|AS.).?', r.type),
    soup.residues())
print "Charged residues:"
for res in residues:
  print res.tag()+'-'+res.type



# find contact residue pairs by decomposing comparisons
# into functions that can be used later

def in_contact(res1, res2):
  for a1 in res1.atoms():
    for a2 in res2.atoms():
      if v3.distance(a1.pos, a2.pos) < 4:
        return True
  return False

print "Contact pairs:"
residues = []
for res in soup.residues():
  if res.type not in data.solvent_res_types:
    residues.append(res)
n = len(residues)
for i in range(n):
  for j in range(i+3, n):
    res1 = residues[i]
    res2 = residues[j]
    if in_contact(res1, res2):
      print (residues[i].tag(), residues[j].tag())



# find residues around catalytic sites
# display the results using:
#    >> pdbshow -b pdb/1be9.binding.pdb

ligand_residues = [r for r in soup.residues() if r.chain_id == "B"]

binding_residues = []
for res in soup.residues():
  if res not in ligand_residues and res.type not in data.solvent_res_types:
    for lig_res in ligand_residues:
      if in_contact(res, lig_res):
        binding_residues.append(res)

for a in soup.atoms():
  a.bfactor = 0.0
for res in binding_residues:
  for a in res.atoms():
    a.bfactor = 1.0

soup.write_pdb('1be9.binding.pdb')



# find hydrogen bonds between N-H and O atoms using a dirty
# method, looking at only heavy atom and distance criteria

def neighbours(atom, test_atoms, distance):
  results = []
  for test_atom in test_atoms:
    if v3.distance(test_atom.pos, atom.pos) < distance:
      if abs(test_atom.res_num - atom.res_num) > 2:
        results.append(test_atom)
  return results

atoms = soup.atoms()
nitrogens = filter(lambda a: a.element=="N", atoms)
oxygens = filter(lambda a: a.element=="O", atoms)

print "Hydrogen bonds:"
for oxygen in oxygens:
  bonded_nitrogens = neighbours(oxygen, nitrogens, 3.5)
  for nitrogen in bonded_nitrogens:
    d = v3.distance(oxygen.pos, nitrogen.pos)
    print (str(oxygen), str(nitrogen)), "%.2f" % d 



# find surface residues/set b-factor using the inbuilt
# accessible surface-area calculator with a water probe radius
# of 1.4 angstroms. writes out PDB file to be viewed:
#     >> pdbshow -b 1be9.asa.pdb
# And:
#     >> pdbshow -b 1be9.res.asa.pdb

atoms = []
for res in soup.residues():
  if res.type not in data.solvent_res_types:
    atoms.extend(res.atoms())
pdbatoms.add_radii(atoms)
areas = asa.calculate_asa(atoms, 1.4)
for a, area in zip(atoms, areas):
  a.asa = area
  a.bfactor = area
for a in soup.atoms():
  if not hasattr(a, 'asa'):
    a.bfactor = 0.0
print "Making atomic ASA pdb", '1be9.asa.pdb'
soup.write_pdb('1be9.asa.pdb')

print "Making residue ASA pdb", '1be9.asa.pdb'
for res in soup.residues():
  sum_b = sum(a.bfactor for a in res.atoms())
  for a in res.atoms():
    a.bfactor = sum_b
soup.write_pdb('1be9.res.asa.pdb')



# find salt bridges using distances between charged N atoms
# of the two positive amino acids, and the O atoms of the
# two negative amino acids

def get_atom_from_res(res, res_type, atom_type):
  if res_type == res.type and res.has_atom(atom_type):
    return res.atom(atom_type)
  return None

nitrogens = []
for res in soup.residues():
  nitrogens.append(get_atom_from_res(res, 'LYS', 'NZ'))
  nitrogens.append(get_atom_from_res(res, 'ARG', 'NH1'))
  nitrogens.append(get_atom_from_res(res, 'ARG', 'NH2'))
nitrogens = filter(lambda a:a is not None, nitrogens)

oxygens = []
for res in soup.residues():
  oxygens.append(get_atom_from_res(res, 'ASP', 'OD1'))
  oxygens.append(get_atom_from_res(res, 'ASP', 'OD2'))
  oxygens.append(get_atom_from_res(res, 'GLU', 'OE1'))
  oxygens.append(get_atom_from_res(res, 'GLU', 'OE2'))
oxygens = filter(lambda a:a is not None, oxygens)

bridges = set()
for nitrogen in nitrogens:
  for oxygen in neighbours(nitrogen, oxygens, 4):
    bridges.add((nitrogen.res_tag(), oxygen.res_tag()))

print "Salt bridges:"
for bridge in bridges:
  print bridge



# Combine the ASA with salt bridge calculations to
# quickly identify which salt bridge is on the surface

print "Salt bridge exposure area:"
for bridge in bridges:
  tag1, tag2 = bridge
  res1 = soup.residue_by_tag(tag1)
  res2 = soup.residue_by_tag(tag2)
  asa1 = sum(a.asa for a in res1.atoms())
  asa2 = sum(a.asa for a in res2.atoms())
  print bridge, asa1 + asa2




