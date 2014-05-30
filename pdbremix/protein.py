# encoding: utf-8

__doc__ = """

Common protein manipulations.

This is a collection of some common PDB and Soup manipulations
of protein strucutres.
"""


import string

import data
import pdbatoms
import util
import pdbtext
import v3


def strip_solvent_pdb(pdb):
  new_pdb = util.fname_variant(pdb)
  txt = pdbtext.strip_solvent(open(pdb).read())
  open(new_pdb, 'w').write(txt)
  return new_pdb


def find_ca_of_resname(atoms, resname):
  for atom in atoms:
    if pdbatoms.split_tag(resname) == \
       (atom.chain_id, atom.res_num, atom.res_insert):
      if "CA" == atom.type:
        return atom
  raise IndexError, "Can't find atom %s" % resname


def get_pdb_transform(pdb, center_res, top_res):
  """
  Returns a transformation matrix that centers pdb to 
  center_res on the z-axis and moves top_res above center_res
  on the y-axis
  """
  soup = pdbatoms.Soup(pdb)
  atoms = soup.atoms()
  soup_center = pdbatoms.get_center(atoms)
  translation = v3.translation(-soup_center)
  soup.transform(translation)
  result = translation

  view = v3.vector(0, 0, 1)

  if center_res is not None:
    center_atom = find_ca_of_resname(soup.atoms(), center_res)
    axis = v3.cross(view, center_atom.pos)
    angle = v3.vec_dihedral(view, axis, center_atom.pos)
    rotation1 = v3.rotation(axis, angle)
    soup.transform(rotation1)
    result = v3.combine(rotation1, result)

  if top_res is not None:
    top_atom = find_ca_of_resname(soup.atoms(), top_res)
    top_dir = v3.vector(0, 1, 0)
    axis = v3.vector(view)
    angle = v3.vec_dihedral(top_dir, axis, top_atom.pos)
    rotation2 = v3.rotation(axis, angle)
    result = v3.combine(rotation2, result)

  return result


def transform_pdbs_to_residues_of_first_pdb(pdbs, center_res, top_res):
  transform = get_pdb_transform(pdbs[0], center_res, top_res)
  new_pdbs = []
  for pdb in pdbs:
    new_pdb = util.fname_variant(pdb)
    soup = pdbatoms.Soup(pdb)
    soup.transform(transform)
    soup.write_pdb(new_pdb)
    new_pdbs.append(new_pdb)
  return new_pdbs


def transformed_soup_from_pdb(
    pdb, center_res=None, top_res=None, 
    width=None, height=None, frame_residues=None):
  soup = pdbatoms.Soup(pdb)
  if center_res and top_res:
    transform = get_pdb_transform(pdb, center_res, top_res)
    soup.transform(transform)
  if frame_residues:
    resnames = [pymol_id_from_res_tag(r) for r in frame_residues]
    soup.frame_pymol_script = "zoom (%s)\n" % ' or '.join(resnames)
  if width: soup.width = width
  if height: soup.height = height
  return soup


def is_peptide_connected(res1, res2):
  if res1.has_atom('CA') and \
     res1.has_atom('C') and \
     res2.has_atom('N') and \
     res2.has_atom('CA'):
     d = v3.distance(res1.atom('C').pos, res2.atom('N').pos)
     if d < 2.0:
       return True
  return False


def find_chains(soup):
  residues = soup.residues()
  n = len(residues)
  if n == 0:
    return
  i_chain = 0
  chain_id = string.ascii_uppercase[i_chain]
  for i in range(0, n):
    res = residues[i]
    if i == 0:
      is_connected_to_prev = False
    else:
      prev_res = residues[i-1]
      is_connected_to_prev = is_peptide_connected(prev_res, res)
    if i < n-1:
      next_res = residues[i+1]
      is_connected_to_next = is_peptide_connected(res, next_res)
    else:
      is_connected_to_next = False
    if is_connected_to_prev or is_connected_to_next:
      res.set_chain_id(chain_id)
      if not is_connected_to_next:
        i_chain += 1
        chain_id = string.ascii_uppercase[i_chain]


