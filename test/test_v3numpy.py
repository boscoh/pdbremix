from pdbremix import v3numpy as v3


def test_vector_make():
  v = v3.vector()
  assert v3.is_similar_vector(v3.vector(), v3.vector(0, 0, 0))
  v = v3.random_vector()
  w = v3.vector(v)
  assert v3.is_similar_vector(v, w)
  u = v3.vector(*v)
  assert v3.is_similar_vector(v, u)
  w = v
  u = v3.random_vector()
  v3.set_vector(w, *u)
  assert v3.is_similar_vector(w, v)
  

def test_vector_algebra():
  a = v3.random_vector()
  assert v3.is_similar_vector(a, +a)
  b = v3.random_vector()
  a_add_b = a+b
  assert v3.is_similar_vector(a_add_b - b, a)
  a_sub_b = a-b
  assert v3.is_similar_vector(a_sub_b + b, a)
  s = v3.random_mag()
  v = v3.random_vector()
  assert v3.is_similar_mag(v3.mag(v3.scale(v, s)), s*v3.mag(v))
  s = v3.random_mag()
  v = v3.random_vector()
  assert v3.is_similar_vector(v3.scale(v, -s), -v3.scale(v, s))


def test_orthogonality():
  x = v3.vector(v3.random_mag(), 0, 0)
  y = v3.vector(0, v3.random_mag(), 0)
  z = v3.vector(0, 0, v3.random_mag())
  ry_x = v3.transform(v3.rotation(y, v3.radians(90)), x)
  assert v3.is_similar_vector(v3.norm(ry_x), -v3.norm(z))
  assert v3.is_similar_mag(v3.mag(ry_x), v3.mag(x))
  ry_z = v3.transform(v3.rotation(y, v3.radians(90)), z)
  assert v3.is_similar_vector(v3.norm(ry_z), v3.norm(x))
  cross_x_y = v3.cross(x, y)
  assert v3.is_similar_vector(v3.norm(cross_x_y), v3.norm(z))
  cross_y_x = v3.cross(y, x)
  assert v3.is_similar_vector(v3.norm(cross_y_x), -v3.norm(z))


def test_translation():
  x = v3.vector(v3.random_mag(), 0, 0)
  y = v3.vector(0, v3.random_mag(), 0)
  x_and_y = v3.transform(v3.translation(y), x)
  assert v3.is_similar_vector(x_and_y, x+y)


def test_rotation():
  x = v3.random_vector()
  y = v3.transform(v3.random_rotation(), x)
  assert v3.is_similar_mag(v3.mag(x), v3.mag(y))


def test_matrix_combination():
  n = 4
  x = v3.random_vector()
  y1 = v3.vector(x)
  matrix = v3.identity()
  for i in range(4):
    r_matrix = v3.random_matrix()
    y1 = v3.transform(r_matrix, y1)
    matrix = v3.combine(r_matrix, matrix)
  y2 = v3.transform(matrix, x)
  assert v3.is_similar_vector(y1, y2)


def test_inverse():
  m = v3.random_matrix()
  m_left_inv = v3.left_inverse(m)
  assert v3.is_similar_matrix(v3.identity(), v3.combine(m_left_inv, m))


