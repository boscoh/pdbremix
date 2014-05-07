import math
import random
from array import array


def radians(degrees):
  return math.pi / 180.0 * degrees


def degrees(radians):
  return 180.0 / math.pi*degrees * radians


class Vec3(array):
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
  return Vec3(s*v[0], s*v[1], s*v[2])


def vector(*args):
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


class Matrix3d:

  def __init__(self):
    self.elem00 = 1.0
    self.elem01 = 0.0
    self.elem02 = 0.0
    self.elem03 = 0.0

    self.elem10 = 0.0
    self.elem11 = 1.0
    self.elem12 = 0.0
    self.elem13 = 0.0

    self.elem20 = 0.0
    self.elem21 = 0.0
    self.elem22 = 1.0
    self.elem23 = 0.0

    self.elem30 = 0.0
    self.elem31 = 0.0
    self.elem32 = 0.0
    self.elem33 = 1.0

  def __str__(self):
    row1 = "  [% .2f, % .2f, % .2f ]\n" % \
              (self.elem00, self.elem01, self.elem02)
    row2 = "  [% .2f, % .2f, % .2f ]\n" % \
              (self.elem10, self.elem11, self.elem12)
    row3 = "  [% .2f, % .2f, % .2f ]\n" % \
              (self.elem20, self.elem21, self.elem22)
    row4 = "  [ ------------------ ]\n"
    row5 = "  [% .2f, % .2f, % .2f ]" % \
              (self.elem30, self.elem31, self.elem32)
    return row1 + row2 + row3 + row4 + row5

  def elem(self, i, j, val=None):
    if val is not None:
      if j==0:
        if i==0: self.elem00 = val
        if i==1: self.elem10 = val
        if i==2: self.elem20 = val
        if i==3: self.elem30 = val
      if j==1:
        if i==0: self.elem01 = val
        if i==1: self.elem11 = val
        if i==2: self.elem21 = val
        if i==3: self.elem31 = val
      if j==2:
        if i==0: self.elem02 = val
        if i==1: self.elem12 = val
        if i==2: self.elem22 = val
        if i==3: self.elem32 = val
      if j==3:
        if i==0: self.elem03 = val
        if i==1: self.elem13 = val
        if i==2: self.elem23 = val
        if i==3: self.elem33 = val

    if j==0:
      if i==0: return self.elem00
      if i==1: return self.elem10
      if i==2: return self.elem20
      if i==3: return self.elem30
    if j==1:
      if i==0: return self.elem01
      if i==1: return self.elem11
      if i==2: return self.elem21
      if i==3: return self.elem31
    if j==2:
      if i==0: return self.elem02
      if i==1: return self.elem12
      if i==2: return self.elem22
      if i==3: return self.elem32
    if j==3:
      if i==0: return self.elem03
      if i==1: return self.elem13
      if i==2: return self.elem23
      if i==3: return self.elem33


def matrix_elem(matrix, i, j, val=None):
  if val is not None:
    if j==0:
      if i==0: matrix.elem00 = val
      if i==1: matrix.elem10 = val
      if i==2: matrix.elem20 = val
      if i==3: matrix.elem30 = val
    if j==1:
      if i==0: matrix.elem01 = val
      if i==1: matrix.elem11 = val
      if i==2: matrix.elem21 = val
      if i==3: matrix.elem31 = val
    if j==2:
      if i==0: matrix.elem02 = val
      if i==1: matrix.elem12 = val
      if i==2: matrix.elem22 = val
      if i==3: matrix.elem32 = val
    if j==3:
      if i==0: matrix.elem03 = val
      if i==1: matrix.elem13 = val
      if i==2: matrix.elem23 = val
      if i==3: matrix.elem33 = val

  if j==0:
    if i==0: return matrix.elem00
    if i==1: return matrix.elem10
    if i==2: return matrix.elem20
    if i==3: return matrix.elem30
  if j==1:
    if i==0: return matrix.elem01
    if i==1: return matrix.elem11
    if i==2: return matrix.elem21
    if i==3: return matrix.elem31
  if j==2:
    if i==0: return matrix.elem02
    if i==1: return matrix.elem12
    if i==2: return matrix.elem22
    if i==3: return matrix.elem32
  if j==3:
    if i==0: return matrix.elem03
    if i==1: return matrix.elem13
    if i==2: return matrix.elem23
    if i==3: return matrix.elem33


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
         val += a.elem(k,i) * b.elem(j,k)
      c.elem(j, i, val)
    # c(3,i) is the translation vector
    val = a.elem(3,i)
    for k in range(0,3):
      val += a.elem(k,i) * b.elem(3,k)
    c.elem(3, i, val)
  return c


def identity():
  return Matrix3d()


def left_inverse(m):
  rot_t = identity()
  for i in range(0, 3):
    for j in range(0, 3):
      val = m.elem(j,i)
      rot_t.elem(i,j,val)
  t = vector(m.elem(3,0), m.elem(3,1), m.elem(3,2))
  t_inv = transform(rot_t, t)
  for i in range(0, 3):
    rot_t.elem(3, i, -t_inv[i])
  return rot_t


def rotation(axis, theta):
  """ matrix to rotate a vector at origin"""
  m = identity()
  c = math.cos(theta)
  s = math.sin(theta)
  t = 1.0 - c
  x, y, z = norm(axis)
  m.elem(0, 0, t*x*x       + c)
  m.elem(0, 1, t*x*y + z*s    )
  m.elem(0, 2, t*x*z - y*s    )
  m.elem(1, 0, t*y*x - z*s    )
  m.elem(1, 1, t*y*y       + c)
  m.elem(1, 2, t*y*z + x*s    )
  m.elem(2, 0, t*z*x + y*s    )
  m.elem(2, 1, t*z*y - x*s    )
  m.elem(2, 2, t*z*z       + c)
  return m


def translation(t):
  """ matrix to translate a vector"""
  m = identity()
  for i in range(3):
    m.elem(3, i, t[i])
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
      if not is_similar_mag(a.elem(i,j), b.elem(i,j), small):
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




