import math
import random

import numpy as np

def is_similar_mag(a, b, small=1E-5):
  return abs(abs(a)-abs(b)) <= small


def vector(*args):
  n_arg = len(args)
  if n_arg == 0:
    return np.zeros(3, dtype=np.float)
  if n_arg == 1:
    data = args[0]
    if len(data) == 3:
      return np.array(data, copy=True)
    raise TypeError('vector() with 1 argument must have 3 elements')
  if n_arg == 3:
    return np.array(args, dtype=np.float, copy=True)
  else:
    raise TypeError('vector() takes 0, 1 or 3 arguments')


def set_vector(*args):
  "Changes values of a vector in place"
  vector = args[0]
  if len(args) == 2:
    vector[:] = args[1]
  elif len(args) == 4:
    vector[:] = args[1:4]


def crd(vector):
  "Returns values of vector as a sequence of floats"
  return vector


def is_similar_matrix(m1, m2, small=1E-5):
  iter1 = np.ndenumerate(m1)
  iter2 = np.ndenumerate(m2)
  for (i1, val1), (i2, val2) in zip(iter1, iter2):
    if not is_similar_mag(val1, val2, small):
      return False
  return True


is_similar_vector = is_similar_matrix


def mag(vector):
  return np.sqrt(np.dot(vector, vector))


def scale(vector, s):
  return  s*vector


def norm(vector):
  return vector/mag(vector)


radians = np.radians

degrees = np.degrees

cross = np.cross

dot = np.dot


def identity():
  m = np.zeros((4, 3))
  m[:3,:3] = np.eye(3)
  return m


def transform(matrix, vector):
  return np.dot(matrix[:3,:3], vector) + matrix[3,:]  


def left_inverse(matrix):
  inverse = identity()
  r = matrix[:3,:3].transpose()
  inverse[:3,:3] = r
  inverse[3,:] = -np.dot(r, matrix[3,:])
  return inverse


# from http://stackoverflow.com/a/6802723
# uses the right hand screw rule for theta
def rotation(axis, theta):
  m = identity()
  a = np.cos(theta/2)
  b, c, d = norm(axis) * np.sin(theta/2)
  m[0] = [a*a+b*b-c*c-d*d, 2*(b*c-a*d),     2*(b*d+a*c)    ]
  m[1] = [2*(b*c+a*d),     a*a+c*c-b*b-d*d, 2*(c*d-a*b)    ]
  m[2] = [2*(b*d-a*c),     2*(c*d+a*b),     a*a+d*d-b*b-c*c]
  return m


def translation(displacement):
  m = identity()
  m[3,:] = displacement
  return m


def combine(m1, m2):
  m3 = identity()
  m3[:3,:3] = np.dot(m1[:3,:3], m2[:3,:3])
  m3[3,:] = np.dot(m1[:3,:3], m2[3,:]) + m1[3,:]
  return m3




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



