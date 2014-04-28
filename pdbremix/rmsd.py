import v3
import pdbatoms

try:
  import numpy as np
except:
  raise Exception("Can't use RMSD routines without NUMPY.")


def rmsd(in_crds1, in_crds2):
  """Returns RMSD between 2 sets of [nx3] np array"""

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
  return np.sqrt(rmsd_sq)


def superposition(in_ref_crds, in_target_crds):
  """
  Returns best transform m between two sets of [nx3] arrays.
  
  Direction: transform(m, ref_crd) => target_crd.
  """
  ref_crds = np.array(in_ref_crds)
  target_crds = np.array(in_target_crds)
  assert(ref_crds.shape[1] == 3)
  assert(ref_crds.shape == target_crds.shape)

  correlation_matrix = np.dot(np.transpose(ref_crds), target_crds)
  v, s, w = np.linalg.svd(correlation_matrix)
  is_reflection = (np.linalg.det(v) * np.linalg.det(w)) < 0.0
  if is_reflection:
    v[-1,:] = -v[-1,:]

  rot33 = np.dot(v, w)
  m = v3.identity()
  m[:3,:3] = rot33.transpose()

  return m


def sum_rmsd(crds1, crds2):
  sum_squared = 0.0
  for crd1, crd2 in zip(crds1, crds2):
    sum_squared += v3.distance(crd1, crd2)**2
  return np.sqrt(sum_squared/float(len(crds1)))
  

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

  if not transform_pdb1:
    return rmsd(crds1, crds2)

  transform_1_to_2 = superposition(crds1, crds2)
  polymer1.transform(transform_1_to_2)

  polymer1.transform(v3.translation(center2))
  polymer2.transform(v3.translation(center2))

  polymer1.write_pdb(transform_pdb1)
  return sum_rmsd(crds1, crds2)


  