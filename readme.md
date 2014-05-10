# PDBREMIX

`pdbremix` is a python library for computational structural biology.

The library provides a light API; it has no external python dependencies, so works in jython and pypy. 

The library can be split up into:

1. standalone tools for PDB structures
2. wrappers around external tools for PDB structures 
3. tools to analyze MD trajectories
3. python interface to analyze and edit PDB structures
4. python interface to run molecular-dynamics simulations
5. python interface to analyse molecular dynamics trajectories


## Standalone PDB structure analysis tools

`pdbremix` provides a bunch of tools to investigate PDB structures that can be used out of the box:

- `pdbfetch` fetches PDB files from the RCSB website
- `pdbheader` displays summary of PDB files
- `pdbseq` displays sequences in a PDB file
- `pdbchain` extracts chains from a PDB file
- `pdbcheck` checks for common defects in a PDB file

These tools implement standard structural biology algorithms using pure Python:

- `pdbvol` calculates the volume
- `pdbasa` calculates the accessible surface-area
- `pdbrmsd` calculates RMSD between structures

- `md2pdb` converts MD restart files into PDB

You can get help for these with the `-h` option on the command-line. 


## Wrappers around external tools

There are many wonderful tools for computational structural biology that have less-than-stellar interfaces. `pdbremix` wraps some of these tools with a friendlier interface and extra functionality. 

To use these tools, first check what binaries that `pdbremix` can find on your system:

	> checkpdbremix

To hand-code the exact location of the binaries, or to run the binaries with exotic flags, edit the config file with the `-o` option such as:

	> vi `checkpdbremix -o`

## Wrappers around external PDB structure tools

These command-line tools provide a better interface to common tools such as PYMOL, MODELLER and THESEUS:

- `pdbshow` wrapper for PYMOL to display PDB with residue centring and b-factor colouring
- `pdboverlay` wraps MAFFT, THESEUS and PYMOL to display homologous proteins
- `pdbinsert` wraps MODELLER to build loops for gaps in a PDB structure


## Trajectory analysis tools

`pdbremix` provides a simplified interface to run molecular-dynamics on 3 standard MD packages: AMBER11+, GROMACS4.5+ and NAMD2.8+. `pdbremix` assumes that associated files for an MD trajectory *has a common base name*. 

In AMBER:

1. topology: md.top
2. coordinates: md.trj 
3. velocities: md.vel.trj

In GROMACS:

1. topology: md.top and associated md.\*.itp files
2. restart coordinates: md.gro
3. coordinates/velocities: md.trr 

In NAMD:

1. topology: md.psf
2. coordinates: md.dcd 
3. velocities: md.vel.dcd

Once named properly, these tools can be used:

- `trajstep` extracts basic trajectory parameters
- `trajvar` calculates energy and RMSD of trajectory

Use these tools to display trajectories: 

- `trajvmd` open trajectories in VMD *recommended*
- `trajchim` open trajectories in CHIMERA
- `trajpym` open trajectories in PYMOL *AMBER only*

These are some package specific tools: 

- `traj2amb` converts NAMD/GROMACS to AMBER trajectories _without_ solvent
- `grotrim` trim GROMACS .trr trajectory files


## Molecular dynamics tools

The library was originally written to do PUFF steered-molecular dynamics simulation that uses a pulsed force application. This method does not require the underlying MD package to support the method and can be carried out by manipulating restart files. The utilities to do this are:

- `pdbrestart`
- `pdbminimize`
- `pdbequil`
- `puff` runs a PUFF simulation
- `puffshow` displays pulling residues in PUFF sim


## Python Interface to manipulate PDB Structures

### Vector geometry library

For a structural biology library, you need a good vector geometry library. Here we provide `v3`:

	from pdbremix import v3

The library switches between a pure python `v3array` version or if you have a numpy, a numpy-dependent version `v3numpy`. 

If you want specifically the python version:

	import pdbremix.v3array as v3

Or the numpy version: 

	import pdbremix.v3numpy as v3

The interface is accessed through functions rather than methods or operator overloading (except for arithmetic vector). This allows the library swapping. Vectors are created by a function:

		v = v3.vector() # the zero vector

	z = v3.vector(1,2,3)

	w = v3.vector(z) # a copy

As vectors are subclassed from arrays, to access components:

	print v[0], v[1], v[2]

Vector functions return by value except for `set_vector`, which changes components in place:

	v = v3.vector(1, 1, 1)
	set_vector(v, 2, 2, 2)

Here are a set of common vector operations, from which most vector operations can be constructed:

	mag(v)
	Returns magnitude of v.
	scale(v, s)
	Returns the vector of v where components multiplied by s.
	dot(v1, v2)
	Returns the dot product of two vectors.
	cross(v1, v2)
	Returns the cross product of two vectors.
	norm(v)
	Returns vector of v where magnitude is normalised to 1.
	parallel(v, axis)
	Returns vector component of v parallel to axis.
	perpendicular(v, axis)
	Returns vector component of v perpendicular to axis.

We also need a representations for affine (rotation+translation) transforms of vectors. In `v3`, transforms are created by these functions:

	identity()
		Returns the identity transform.
	rotation(axis, theta)
		Returns transform that rotates around the origin.
	rotation_at_center(axis, theta, center)
		Returns transform that rotates around center.
	translation(t)
		Returns transform that translates a vector by t.
	left_inverse(m)
		Returns the left inverse of m.
			Example:
			  combine(left_inverse(m), m) == identity().

Which can be combined into single transforms:

	combine(a, b)
	Returns transform that combines two transforms.

The purpose of transforms is to transform a vector, and this is done by:

	transform(matrix, v)
			Returns transform that translates by displacement vector t.

You typically won’t want to access the elements of a matrix, but if you do, it is through this in-place function:

	matrix_elem(matrix, i, j, val=None)
	Reads/writes the elements of an affine transform.
	1. 3x3 rotational component;
	    matrix_elem(m, i, j) for i=0..3, j=0..3
	2. 3x1 translational component:
	    matrix_elem(m, 3, i) for i=0..3)

distance(p1, p2)
Returns distance between two points
degrees(radians)
Converts radians to degrees, better for reporting.
radians(degrees)
Converts degrees to radians, which is used in math functions.
normalize_angle(angle)
Returns angle in radians that is [-pi, pi]()
vec_angle(a, b)
Returns angle in radians between a and b.
vec_dihedral(a, axis, c)
Returns dihedral angle between a and c, along the axis.
dihedral(p1, p2, p3, p4)
Returns dihedral angle defined by the four positions.

is_similar_mag(a, b, small=0.0001)
Evaluates similar magnitudes to within small.
is_similar_matrix(a, b, small=0.0001)
Evaluates similar matrixes through matrix components.
is_similar_vector(a, b, small=0.0001)
Evaluates similar matrixes through matrix components.

get_center(crds)
Returns the geometric center of a bunch of positions.
get_width(crds)
Returns the maximum width between any two crds in the group.

random_mag()
Returns a random positive number from [0, 90]() for testing.
random_matrix()
Returns a random transformation matrix for testing.
random_real()
Returns a random real +/- from [-90, 90]() for testing.
random_rotation()
Returns a random rotational matrix for testing.
random_vector()
Returns a random vector for testing.

v3.py
v3array.py
v3numpy.py

### Reading in a Soup


The main object for manipulating PDB structures is the Soup object in `pdbatoms`. It is essentially a collection of atoms with a coordinated list of residues. A residue in this case represents a collection of atoms that forms a recognizable chemical group, whether a distinct molecule in the case of solvent, ions and ligands, or an actual residue that forms a polymer, as in amino acids in a protein or a nucleic acid in DNA. 

Specifically chains are separate operations on a Soup and not explicitly represented. 

	 from pdbremix import fetch
	 from pdbremix import pdbatoms
	
	 fetch.fetch_pdb("1be9")
	 soup = pdbatoms.Soup("1be9.pdb")
	
	 print len(soup)
	
	 print soup.residues()
	
	 print soup.atoms()

util.py
asa.py
data.py
fetch.py
pdbatoms.py
pdbtext.py
protein.py
pymol.py
rmsd.py

volume.py

simulate.py
force.py
amber.py
gromacs.py
namd.py

traj.py


## A vector library

- numerics
- function based, not object based
- allows wrapping around other libraries
- not operator overloading, with matrices, it's too easy to get confused
- allows non-numpy library, so pypy
	pypy `which pdbasa` 1be9.pdb

# representing PDB structures

- I have written quite a few variations of the pdbatoms library, I've found the right level of abstraction, the simplest required to do all the things that I find useful
- the PDB really has 
	1. atoms
	2. residues
- chains are really not implemented coherently and should be considered on an ad-hoc basis
- really important to have a PDB structure with a general geometrical transformation library
- move pieces around, splice things around, use proper vector geometry

# patching PDB structures
- delete extra waters
- delete non-standard amino acids
- look for steric clashes
- identify chain breaks/missing amino acids
- reduce
- patch with modeller
- remove extra models
- alternate conformations

# abstracting MD simulations

- coercing all necessary files under the same basename
- making topology, coordinate, velocity files
- choosing a reasonable robust subset of simulation parameters
- can handle proteins
	- only standard amino acids
	- checks for disulfide bonds
	- handle multiple chains
	- ignores cystallographic waters
	- ignores hydrogen atoms
- explicit waters
	- periodic cubic box
	- 10 Å padding
	- Langevin thermometer
	- Nose-Hoover barometer
	- Particle Ewald-Mesh Electrostatics
	- no bond constraints on protein
- implicit solvent
	- Generalized Born electrostatics
	- Surface Area tension hydrophobic term

# PDB fetch

# Surface Area Calculation

# Volume Calculation

# RMSD
- raw calculation, in the only place that numpy is required, it uses the classic SVD decomposition in numpy to calculate the optimal superposition between two sets of points
- pdboverlay

# Interacting with viewers
- frustrating to get the view you want, loading trajectories is a multiple step processing in viewers
- viewers Pymol, Chimera, VMD have scripting languages
- if we enforce our naming convention, can really save time
- e.g. load trajectories using a simple command line
- pipe in useful information and *transformations* 
	`pdbpym -b -c B:8 -t B:7 1be9.pdb`

# Restart files for maximum flexibility

# PUFF steered molecular dynamics

# Reading trajectories

Combines trajectory reading with topology
You need atomic masses and charges for calculations
Reads as coordinates and or as pdbatoms.Polymer structure
Common interface for AMBER, NAMD and GROMACS

# Notes
https://github.com/synapticarbors/pyqcprot
mpiexec -np 40 /home/bosco/bin/gromacs-4.0.7/bin/mdrun -v -s md.tpr -cpi md.cpt -append -deffnm md \>& md.mdrun.restart.log
gromacs ligand http://www.dddc.ac.cn/embo04/practicals/9\_15.htm
installing amber on mac http://amberonmac.blogspot.com.au/
amber ff on gromacs http://www.somewhereville.com/?p=114
ffamber http://ffamber.cnsm.csulb.edu/
http://web.mit.edu/vmd\_v1.9.1/namd-tutorial-unix.pdf
http://bionano.physics.illinois.edu/Tutorials/ssbTutorial.pdf

