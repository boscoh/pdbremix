# encoding: utf-8

__doc__ = """ 

Calculate the accessible surface area of a bunch of atoms

It uses the simple Shrake-Rupley algorithm, that generates a relatively
uniform density of dots over every atoms and eliminates those within the
sphere of another atom. The  remaining dots is used to calculate the area.

Reference: A. Shrake & J. A. Rupley. "Environment and Exposure to Solvent of
Protein Atoms. Lysozyme and Insulin." J Mol Biol. 79 (1973) 351- 371. 
"""

import math

import v3
import pdbatoms


def generate_sphere_points(n):
  """
  Returns list of 3d coordinates of points on a sphere using the
  Golden-Section Spiral algorithm.
  """
  points = []
  inc = math.pi * (3 - math.sqrt(5))
  offset = 2 / float(n)
  for k in range(int(n)):
    y = k * offset - 1 + (offset / 2)
    r = math.sqrt(1 - y*y)
    phi = k * inc
    points.append(v3.vector(math.cos(phi)*r, y, math.sin(phi)*r))
  return points


def find_neighbor_indices(atoms, probe, k):
  """
  Returns list of indices of atoms within probe distance to atom k. 
  """
  neighbor_indices = []
  atom_k = atoms[k]
  radius = atom_k.radius + probe + probe
  indices = range(k)
  indices.extend(range(k+1, len(atoms)))
  for i in indices:
    atom_i = atoms[i]
    dist = v3.distance(atom_k.pos, atom_i.pos)
    if dist < radius + atom_i.radius:
      neighbor_indices.append(i)
  return neighbor_indices


def calculate_asa(atoms, probe, n_sphere_point=960):
  """
  Returns list of accessible surface areas of the atoms,
  using the probe and atom radius to define the surface.
  """
  sphere_points = generate_sphere_points(n_sphere_point)
 
  const = 4.0 * math.pi / len(sphere_points)
  areas = []
  for i, atom_i in enumerate(atoms):
    
    neighbor_indices = find_neighbor_indices(atoms, probe, i)
    n_neighbor = len(neighbor_indices)
    j_closest_neighbor = 0
    radius = probe + atom_i.radius
    
    n_accessible_point = 0
    for point in sphere_points:
      is_accessible = True
      test_point = v3.scale(point, radius) + atom_i.pos
      cycled_indices = range(j_closest_neighbor, n_neighbor)
      cycled_indices.extend(range(j_closest_neighbor))
      
      for j in cycled_indices:
        atom_j = atoms[neighbor_indices[j]]
        r = atom_j.radius + probe
        diff = v3.distance(atom_j.pos, test_point)
        if diff*diff < r*r:
          j_closest_neighbor = j
          is_accessible = False
          break
      if is_accessible:
        n_accessible_point += 1
    
    area = const*n_accessible_point*radius*radius 
    areas.append(area)

  return areas

