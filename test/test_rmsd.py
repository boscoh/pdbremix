from pdbremix import v3
from pdbremix import rmsd
from pdbremix import pdbatoms


def test_superposition():
  n = 4
  crds1 = [v3.random_vector() for i in range(n)]
  random_m = v3.random_rotation()
  crds2 = [v3.transform(random_m, c) for c in crds1]
  rms, m = rmsd.calc_rmsd_rot(crds1, crds2)
  for crd1, crd2 in zip(crds1, crds2):
    assert v3.is_similar_vector(v3.transform(m, crd1), crd2)


def test_rmsd():
  n = 4
  crds1 = [v3.random_vector() for i in range(n)]
  random_m = v3.random_rotation()
  crds2 = [v3.transform(random_m, c) for c in crds1]
  rms, m = rmsd.calc_rmsd_rot(crds1, crds2)
  assert v3.is_similar_mag(0, rms)
  # test the two algorithms separately
  rms, m = rmsd.pyqcprot_rmsd_rot(crds1, crds2)
  assert v3.is_similar_mag(0, rms)
  if rmsd.is_numpy:
    rms, m = rmsd.numpy_svd_rmsd_rot(crds1, crds2)
    assert v3.is_similar_mag(0, rms)


def test_rmsd_pdbatoms():
  p = pdbatoms.Soup('pdb/1ubq.pdb')
  p.transform(v3.random_matrix())
  p.write_pdb('pdb/1ubq-r.pdb')
  rms = rmsd.rmsd_of_pdbs(
    'pdb/1ubq.pdb', 
    'pdb/1ubq-r.pdb', 
    transform_pdb1='pdb/1ubq-q.pdb')
  assert v3.is_similar_mag(0, rms, 1E-2)
  q = pdbatoms.Soup('pdb/1ubq-q.pdb')
  r = pdbatoms.Soup('pdb/1ubq-r.pdb')
  for atom_q, atom_r in zip(q.atoms(), r.atoms()):
    assert v3.is_similar_vector(atom_q.pos, atom_r.pos, 1E-2)


if __name__ == "__main__":
  test_superposition()
  test_rmsd()
  test_rmsd_pdbatoms()