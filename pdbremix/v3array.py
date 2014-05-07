
__doc__ = """
A 3D vector geometry library that works in pure Python.

The vector object is a subclass of array.array. For optimal 
flexibility, the interface is all done through functions, and
not methods. This way we can easily swap out other internal
representations of vectors and matrices, such as in v3numpy.
"""

import math
import random
from array import array


def radians(degrees):
  """
  Converts degrees to radians, which is used in math functions.
  """
  return math.pi / 180.0 * degrees


def degrees(radians):
  """
  Converts radians to degrees, which is better for reporting.
  """
  return 180.0 / math.pi*degrees * radians


class Vec3(array):
  """
  This is a subclass of array to emulate a 3D vector.

  Arithmetic operators overloading are defined to 
  carry out basic artihmetic vector operations. Scaling
  and other non-arithmetic operators are carried out
  with functions.
  """
  def __new__(cls, x=0, y=0, z=0):
    return array.__new__(cls, 'd', (x,y,z))

  def __repr__(self):
    return "Vec3(%f, %f, %f)" % (self[0], self[1], self[2])

  def __add__(self, rhs):
    return Vec3(self[0]+rhs[0], self[1]+rhs[1], self[2]+rhs[2])

  def __iadd__(self, rhs):
    self[0] += rhs[0]
    self[1] += rhs[1]
    self[2] += rhs[2]
    return self

  def __isub__(self, rhs):
    self[0] -= rhs[0]
    self[1] -= rhs[1]
    self[2] -= rhs[2]
    return self

  def __sub__(self, rhs):
    return Vec3(self[0]-rhs[0], self[1]-rhs[1], self[2]-rhs[2])

  def __neg__(self):
    return Vec3(-self[0], -self[1], -self[2])

  def __pos__(self):
    return self


def scale(v, s):
  """
  Returns a vector that has been scaled by s.
  """
  return Vec3(s*v[0], s*v[1], s*v[2])


def vector(*args):
  """
  Returns a new vector, whether zero, or copied from args.
  """
  if len(args) == 0:
    return Vec3(0, 0, 0)
  elif len(args) == 1:
    v = args[0]
    return Vec3(v[0], v[1], v[2])
  elif len(args) == 3:
    return Vec3(*args)
  else:
    raise TypeError('vector() can take 0, 1, 3 arguments only')


def set_vector(*args):
  "Changes values of a vector in place"
  v = args[0]
  if len(args) == 2:
    w = args[1]
    v[:] = w
  elif len(args) == 4:
    v[0], v[1], v[2] = args[1:4]


def mag(v):
  x, y, z = v
  return math.sqrt(x*x + y*y + z*z)


def dot(v1, v2):
  x1, y1, z1 = v1
  x2, y2, z2 = v2
  return x1*x2 + y1*y2 + z1*z2


def norm(v):
  return scale(v, 1.0/mag(v))


def cross(v1, v2):
  x1, y1, z1 = v1
  x2, y2, z2 = v2
  return vector(
      y1*z2 - z1*y2,
      z1*x2 - x1*z2,
      x1*y2 - y1*x2)


class Matrix3d(array):
  """
  Class to represent affine transforms in 3D space.

  This requires a 3x3 matrix and a translation vector.
  Here we choose matrix(3,i) as the translation vector.
  We represent Matrix3d as an array of 12 floats, to be
  accessed by matrix_elem().
  """
  def __new__(cls, matrix_array=(1,0,0,0,1,0,0,0,1,0,0,0)):
    """
    Initializes to the identity matrix.
    """
    return array.__new__(cls, 'd', matrix_array)

  def __str__(self):
    def str3(x, y, z): 
      return "  [% .2f, % .2f, % .2f ]\n" % (x, y, z)
    s = ""
    s += str3(*self[:3])
    s += str3(*self[3:6])
    s += str3(*self[6:9])
    s += "  [ ------------------ ]\n"
    s += str3(*self[9:12])
    return s


def matrix_elem(matrix, i, j, val=None):
  k = i*3 + j
  if val is None:
    return matrix[k]
  else:
    matrix[k] = val
    return val


def identity():
  return Matrix3d()


def transform(matrix, v):
  v_x, v_y, v_z = v
  x = matrix_elem(matrix, 0, 0) * v_x + \
      matrix_elem(matrix, 1, 0) * v_y + \
      matrix_elem(matrix, 2, 0) * v_z + \
      matrix_elem(matrix, 3, 0)
  y = matrix_elem(matrix, 0, 1) * v_x + \
      matrix_elem(matrix, 1, 1) * v_y + \
      matrix_elem(matrix, 2, 1) * v_z + \
      matrix_elem(matrix, 3, 1)
  z = matrix_elem(matrix, 0, 2) * v_x + \
      matrix_elem(matrix, 1, 2) * v_y + \
      matrix_elem(matrix, 2, 2) * v_z + \
      matrix_elem(matrix, 3, 2)
  return vector(x, y, z)


def combine(a, b):
  c = identity()
  for i in range(0, 3):
    for j in range(0, 3):
      val = 0.0
      for k in range(0, 3):
         val += matrix_elem(a,k,i) * matrix_elem(b,j,k)
      matrix_elem(c, j, i, val)
    val = matrix_elem(a,3,i)
    for k in range(0,3):
      val += matrix_elem(a,k,i) * matrix_elem(b,3,k)
    matrix_elem(c, 3, i, val)
  return c


def left_inverse(m):
  rot_t = identity()
  for i in range(0, 3):
    for j in range(0, 3):
      val = matrix_elem(m,j,i)
      matrix_elem(rot_t,i,j,val)
  t = vector(matrix_elem(m,3,0), matrix_elem(m,3,1), matrix_elem(m,3,2))
  t_inv = transform(rot_t, t)
  for i in range(0, 3):
    matrix_elem(rot_t, 3, i, -t_inv[i])
  return rot_t


def rotation(axis, theta):
  """ matrix to rotate a vector at origin"""
  m = identity()
  c = math.cos(theta)
  s = math.sin(theta)
  t = 1.0 - c
  x, y, z = norm(axis)
  matrix_elem(m, 0, 0, t*x*x       + c)
  matrix_elem(m, 0, 1, t*x*y + z*s    )
  matrix_elem(m, 0, 2, t*x*z - y*s    )
  matrix_elem(m, 1, 0, t*y*x - z*s    )
  matrix_elem(m, 1, 1, t*y*y       + c)
  matrix_elem(m, 1, 2, t*y*z + x*s    )
  matrix_elem(m, 2, 0, t*z*x + y*s    )
  matrix_elem(m, 2, 1, t*z*y - x*s    )
  matrix_elem(m, 2, 2, t*z*z       + c)
  return m


def translation(t):
  """ matrix to translate a vector"""
  m = identity()
  for i in range(3):
    matrix_elem(m, 3, i, t[i])
  return m


def is_similar_mag(a, b, small=1E-4):
  return abs(a-b) <= small


def is_similar_vector(a, b, small=1E-4):
  for a_x, b_x in zip(a, b):
    if not is_similar_mag(a_x, b_x, small):
      return False
  return True


def is_similar_matrix(a, b, small=1E-4):
  for i in range(0, 3):
    for j in range(0, 3):
      if not is_similar_mag(matrix_elem(a,i,j), matrix_elem(b,i,j), small):
        return False
  return True


# common to v3list and v3numpy


def parallel(v, axis):
  l = mag(axis)
  if is_similar_mag(l, 0):
    return v
  else:
    return scale(axis, dot(v, axis)/l/l) 


def perpendicular(v, axis):
  return v - parallel(v, axis)


def normalize_angle(angle):
  while abs(angle) > math.pi:
    if angle > math.pi:
      angle -= math.pi*2
    if angle < -math.pi:
      angle += 2*math.pi
  if is_similar_mag(abs(angle + math.pi), 0):
    angle = math.pi
  return angle


def vec_angle(a, b):
  a_len = mag(a)
  b_len = mag(b)
  if is_similar_mag(a_len, 0) or is_similar_mag(b_len, 0):
    return 0.0
  c = dot(a, b) / a_len / b_len
  if c >=  1.0:
    return 0.0
  elif c <= -1.0:
    return math.pi
  else:
    return math.acos(c)  


def vec_dihedral(a, axis, c):
  ap = perpendicular(a, axis)
  cp = perpendicular(c, axis)
  angle = vec_angle(ap, cp)
  if dot(cross(ap, cp), axis) > 0:
    angle = -angle
  return angle


def dihedral(p1, p2, p3, p4):
  return vec_dihedral(p1-p2, p2-p3, p4-p3)


def distance(p1, p2):
  return mag(p1 - p2)


def rotation_at_center(axis, theta, center):
  t = translation(-center)
  r = rotation(axis, theta)
  t_inv = translation(center)
  return combine(t_inv, combine(r, t))


def get_center(crds):
  center = vector()
  for crd in crds:
    center += crd
  return scale(center, 1.0/float(len(crds)))


def get_width(crds):
  center = get_center(crds)
  max_diff = 0
  for crd in crds:
    diff = v3.distance(crd, center)
    if diff > max_diff:
      max_diff = diff
  return 2*max_diff


def random_mag(): 
  return random.uniform(0, 90)


def random_real(): 
  return random.uniform(-90, 90)


def random_vector():
  return vector(random_real(), random_real(), random_real())


def random_rotation():
  return rotation(random_vector(), radians(random_real()))


def random_matrix():
  return combine(random_rotation(), translation(random_vector()))




