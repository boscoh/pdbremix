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



# setup

util.goto_dir('pdb')
fetch.get_pdbs_with_http('1be9', '1ry4')
soup = pdbatoms.Soup('1be9.pdb')



# select polar/charged residues

print "Charged residues:"
charged_residues = []
for res in soup.residues():
  if re.match('(LY.|AR.|GL.|AS.)', res.type[:3]):
    charged_residues.append(res)
for res in charged_residues:
  print res.tag()+'-'+res.type



# contact residues

def in_contact(res1, res2):
  for a1 in res1.atoms():
    for a2 in res2.atoms():
      if v3.distance(a1.pos, a2.pos) < 4:
        return True
  return False

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
      print "Contact", residues[i].tag(), residues[j].tag()



# find residues around catalytic sites

ligand_residues = []
for res in soup.residues(): 
  if res.chain_id == "B":
    ligand_residues.append(res)

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



# find hydrogen bonds

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



# find salt bridges

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



# find surface residues/set b-factor

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



# exposure of salt bridges

print "Salt bridge exposure area:"
for bridge in bridges:
  tag1, tag2 = bridge
  res1 = soup.residue_by_tag(tag1)
  res2 = soup.residue_by_tag(tag2)
  asa1 = sum(a.asa for a in res1.atoms())
  asa2 = sum(a.asa for a in res2.atoms())
  print bridge, asa1 + asa2






