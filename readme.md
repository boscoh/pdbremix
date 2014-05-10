# PDBREMIX

`pdbremix` is a python library for computational structural biology.

The library provides a light API.

; it has no external python dependencies, so works in jython and pypy. 

The library can be split up into:

1. standalone tools for analysis of PDB structures
2. wrappers around external tools for PDB structures 
3. tools to analyze MD trajectories
3. python interface to analyze and edit PDB structures
4. python interface to run molecular-dynamics simulations
5. python interface to analyse trajectories


## Standalone PDB structure analysis tools

`pdbremix` provides a bunch of tools to investigate PDB structures that can be used out of the box:

- `pdbfetch` fetches PDB files from the RCSB website
- `pdbheader` displays summary of PDB files
- `pdbseq` displays sequences in a PDB file

- `pdbcheck` checks for structural defects in a PDB file
- `pdbstrip` to clean up a PDB file
- `pdbchain` extracts chains from a PDB file

The following tools implement standard structural biology algorithms using pure Python:

- `pdbvol` calculates the volume
- `pdbasa` calculates the accessible surface-area
- `pdbrmsd` calculates RMSD between structures

For these tools, you can get a large speed gain if you run them  through `pypy`, where as an example:

	> pypy `which pdbvol` 1be9.pdb

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


## Vector geometry library

As in any structural biology library, we provide a vector geometry library, which is called `v3`:

	from pdbremix import v3

`v3` was designed to be function-based, with vector and transform objects subclassed from arrays. This allows the library to easily switch between a pure Python version and a numpy-dependent version. If you want specifically the python version:

	import pdbremix.v3array as v3

Or the numpy version: 

	import pdbremix.v3numpy as v3

Vectors are created and copied by the `vector` function:

	v = v3.vector() # the zero vector
	z = v3.vector(1,2,3)
	w = v3.vector(z) # a copy

As vectors are subclassed from arrays, to access components:

	print v[0], v[1], v[2]

Most functions return by value, except for `set_vector`, which changes components in place:

	v3.set_vector(v, 2, 2, 2)

Here are a set of common vector operations:

	mag(v)
		Returns magnitude of v.
	scale(v, s)
		Returns the vector of v where components multiplied by s.
	dot(v1, v2)
		Returns the dot product of two vectors.
	cross(v1, v2)
		Returns the cross product of two vectors.
	norm(v)
		Returns vector of v where magnitude is normalised to 1.
	parallel(v, axis)
		Returns vector component of v parallel to axis.
	perpendicular(v, axis)
		Returns vector component of v perpendicular to axis.

Vectors will be used to represent coordinates/points, velocities, displacements etc. Functions are provided to measure their geometric properties:

	distance(p1, p2)
		Returns distance between two points

	vec_angle(a, b)
		Returns angle in radians between a and b.
	vec_dihedral(a, axis, c)
		Returns dihedral angle between a and c, along the axis.
	dihedral(p1, p2, p3, p4)
		Returns dihedral angle defined by the four positions.

	normalize_angle(angle)
		Returns angle in radians that is [-pi, pi]()
	degrees(radians)
		Converts radians to degrees, better for reporting.
	radians(degrees)
		Converts degrees to radians, which is used in math functions.

	get_center(crds)
		Returns the geometric center of a bunch of positions.
	get_width(crds)
		Returns the maximum width between any two crds in the group.

#### Affine Transforms

We also need a representations for affine transforms, which involve a rotation and a translation. 

The purpose of a transform `matrix` is to transform a vector `v`:

	v3.transform(matrix, v)

In `v3`, we represent the transform as a 4x3 matrix with two parts:

1. 3x3 rotational component:

		matrix_elem(m, i, j) for i=0..3, j=0..3

2. 3x1 translational component:

		matrix_elem(m, 3, i) for i=0..3

To read/write the elements of a transform, we use:

	matrix_elem(matrix, i, j, val=None)

Most of the time, you would build a transform from these basic transform generating functions:

	v3.identity()
	v3.rotation(axis, theta)
	v3.rotation_at_center(axis, theta, center)
	v3.translation(t)
	v3.left_inverse(m)

And combine them in the correct sequence:

	c = v3.combine(a, b)

#### Testing Vectors and Transforms

Finally, we introduce functions to test similarity for vectors and transforms:
	
	is_similar_mag(a, b, small=0.0001)
		Evaluates similar magnitudes to within small.
	is_similar_matrix(a, b, small=0.0001)
		Evaluates similar matrixes through matrix components.
	is_similar_vector(a, b, small=0.0001)
		Evaluates similar matrixes through matrix components.

And a bunch of random geometric object generators, used for testing:	

	random_mag()
		Returns a random positive number from [0, 90]() for testing.
	random_matrix()
		Returns a random transformation matrix for testing.
	random_real()
		Returns a random real +/- from [-90, 90]() for testing.
	random_rotation()
		Returns a random rotational matrix for testing.
	random_vector()
		Returns a random vector for testing.

### Reading in a Soup

Here, we look at how to manipulate PDB structure. First, let's grab a PDB structure from the website using the `fetch` module:

	from pdbremix import fetch
	fetch.get_pdbs_with_http(['1be9'])

The main object for manipulating PDB structures is the Soup object in `pdbatoms`. We can read a Soup from a PDB file:

	from pdbremix import pdbatoms
	soup = pdbatoms.Soup("1be9.pdb")

**List of Atoms.** A Soup is essentially a collection of atoms, which we can grab by:

	atoms = soup.atoms()

An atom has attributes:

  - pos (v3.vector)
  - vel (v3.vector)
  - mass (float)
  - charge (float)
  - type (str)
  - element (str)
  - num (int)
  - chain_id (str)
  - res_type (str)
  - res_num (str)
  - res_insert (str)
  - bfactor (float)
  - occupancy (float)
  - alt_conform (str)
  - is_hetatm (bool)

Since `atom.pos` and `atom.vel` are vectors, these can all be manipulated by functions from the `v3` library. Atom contains one special method `transform` that handle affine transforms. 

For example, to move an atom 1.0 angstrom along the X-axis:

	t = v3.translation(v3.vector(1,0,0))
	atom = atoms[0]
	atom.transform(t)

**List of Residues** As well, Soup contains a list of residues:

	residues = soup.residues()

A residue in this case represents a collection of atoms that forms a recognizable chemical group, whether a distinct molecule in the case of solvent, ions and ligands, or an actual residue that forms a polymer, as in amino acids in a protein or a nucleic acid in DNA. 

One useful heuristic is that 
Although chains are representedSpecifically chains are separate operations on a Soup and not explicitly represented. 

	print len(soup)
	nnprint soup.residues()
	
	 print soup.atoms()


util.py
asa.py
data.py
fetch.py
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

