import v3
import pdbatoms
import math

try:
  import numpy as np
  is_numpy = True
except:
  import lib.pyqcprot
  is_numpy = False


def numpy_svd_rmsd_rot(in_crds1, in_crds2):
  """
  Returns best transform m between two sets of [nx3] arrays.
  
  Requires numpy

  Returns RMSD between 2 sets of [nx3] np array
  Direction: transform(m, ref_crd) => target_crd.
  """

  crds1 = np.array(in_crds1)
  crds2 = np.array(in_crds2)
  assert(crds1.shape[1] == 3)
  assert(crds1.shape == crds2.shape)

  n_vec = np.shape(crds1)[0]
  correlation_matrix = np.dot(np.transpose(crds1), crds2)
  v, s, w = np.linalg.svd(correlation_matrix)
  is_reflection = (np.linalg.det(v) * np.linalg.det(w)) < 0.0

  if is_reflection:
    s[-1] = - s[-1]
  E0 = sum(sum(crds1 * crds1)) + \
       sum(sum(crds2 * crds2))
  rmsd_sq = (E0 - 2.0*sum(s)) / float(n_vec)
  rmsd_sq = max([rmsd_sq, 0.0])
  rmsd = np.sqrt(rmsd_sq)

  if is_reflection:
    v[-1,:] = -v[-1,:]
  rot33 = np.dot(v, w)
  m = v3.identity()
  m[:3,:3] = rot33.transpose()

  return rmsd, m


def pyqcprot_rmsd_rot(crds1, crds2):
  rms, rot9 = lib.pyqcprot.calc_rms_rot(crds1, crds2)
  matrix = v3.identity()
  for i in range(3):
    for j in range(3):
       v3.matrix_elem(matrix, i, j, rot9[i*3+j])
  return rms, matrix


def calc_rmsd_rot(crds1, crds2):
  if is_numpy:
    return numpy_svd_rmsd_rot(crds1, crds2)
  else:
    return pyqcprot_rmsd_rot(crds1, crds2)


def sum_rmsd(crds1, crds2):
  sum_squared = 0.0
  for crd1, crd2 in zip(crds1, crds2):
    sum_squared += v3.distance(crd1, crd2)**2
  return math.sqrt(sum_squared/float(len(crds1)))
  

def get_superposable_atoms(polymer, segments, atom_types):
  result = []
  allowed_i = []
  residues = polymer.residues()
  if segments:
    for res_tag_i, res_tag_j in segments:
      i = polymer.get_i_residue(str(res_tag_i))
      j = polymer.get_i_residue(str(res_tag_j))
      if i > j:
        i, j = j, i
      allowed_i.extend(range(i,j+1))
  else:
    allowed_i = range(len(residues))
  for i, residue in enumerate(residues):
    if i in allowed_i:
      result.extend([a for a in residue.atoms()
                     if a.type in atom_types])
  return result


def rmsd_of_pdbs(
    pdb1, pdb2, segments1=[], segments2=[], 
    atom_types=['CA'], transform_pdb1=None):
  polymer1 = pdbatoms.Polymer(pdb1)
  polymer2 = pdbatoms.Polymer(pdb2)

  atoms1 = get_superposable_atoms(polymer1, segments1, atom_types)
  atoms2 = get_superposable_atoms(polymer2, segments2, atom_types)

  crds1 = [a.pos for a in atoms1]
  crds2 = [a.pos for a in atoms2]

  center1 = v3.get_center(crds1)
  center2 = v3.get_center(crds2)

  polymer1.transform(v3.translation(-center1))
  polymer2.transform(v3.translation(-center2))

  rmsd, transform_1_to_2 = calc_rmsd_rot(crds1, crds2)

  if not transform_pdb1:
    return rmsd

  polymer1.transform(transform_1_to_2)

  polymer1.transform(v3.translation(center2))
  polymer2.transform(v3.translation(center2))

  polymer1.write_pdb(transform_pdb1)
  return sum_rmsd(crds1, crds2)


  