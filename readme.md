title: pdbremix documentation
---
---

# pdbremix

`pdbremix` is a library to analyze protein structures and protein simulations

The library consists of:

1. tools to analyze and view PDB structures
2. tools to run MD simulations and analyze MD trajectories
3. python interface to analyze PDB structures
4. python interface for MD simulations and MD trajectories

The library works with PyPy for significant speed-ups.

## Installation

Download from github: 

&nbsp; &nbsp;  [\[zip-package\]](https://github.com/boscoh/pdbremix/archive/master.zip)

And install:

	> python setup.py install

From here, you can access unit tests and example files.

There are many wonderful tools in structural biology that have less-than-stellar interfaces. `pdbremix` wraps these tools to make them easier to use.

To check which tools can be accessed from the path:

	> checkpdbremix

Use the `-o` flag to get the config file for binaries to override (perhaps with exotic flags):

	> vi `checkpdbremix -o`


## Tools to analyze PDB structures

`pdbremix` is a library to analyze PDB structures and MD trajectories. As such, it provides a platform to build command-line tools for PDB files as well as to carry out useful pre-processing of PDB files for external tools.

### Tools in Pure Python

Some of the tools can be used straight out of the box:

- `pdbfetch` fetches PDB files from the RCSB website
- `pdbheader` displays summary of PDB files
- `pdbseq` displays sequences in a PDB
- `pdbchain` extracts chains from a PDB
- `pdbcheck` checks for common defects in a PDB
- `pdbstrip` cleans up PDB for MD simulations

The following tools implement standard algorithms:

- `pdbvol` calculates volume of a PDB
- `pdbasa` calculates accessible surface-area of a PDB
- `pdbrmsd` calculates RMSD between PDB files

For these tools, you get an impressive speed-up if use use `pypy`:

	> pdbfetch 1be9
	> pdbstrip 1be9.pdb
	> pypy `which pdbvol` 1be9.pdb

As usual, detailed help is available with the `-h` flag. 


### Wrappers around External Tools

These following tools wrap external tools to solve some very common (and painful) use-cases in PDB analysis.

- `pdbshow` displays PDB structures in PYMOL with extras.

	PYMOL is a powerful viewer, but it's defaults leave a little to be desired. `pdbshow` runs PYMOL with useful defaults and added functionality:

	  - By default, shows colored chains, ribbons, and sidechains as sticks. 
	  - Define initial viewing frame by a center-residue and a top-residue. Structure is rotated to place the center-residue above the center-of-mass in the middle, and the top-residue above the center-residue.
	  - Color by B-factor using a red-white scale, with limits.
	- Worm mode to show B-factor by variable width
	  - Solvent molecules can be removed, specifically for MD frames that contain too many waters, which will choke PYMOL.  

- `pdboverlay` display homologous PDB files using MAFFT, THESEUS and PYMOL.

	One of the most beautiful results of structural biology is the structural alignment of homologous proteins. `pdboverlay` performs this complex process in one easy step starting from PDB structures:

	- Write fasta sequences from PDB.
	- Align sequences with MAFTT to find homologous regions.
	- Structurally align homologous regions with THESEUS.
	- Display structurally-aligned PDBs using special PYMOL script.  

- `pdbinsert` fill gaps in PDB with MODELLER

	Gaps in PDB structures cause terrible problems in MD simulations. The standard tool to patch gaps is MODELLER, which requires a ton of boilerplate. `pdbinsert` does all the dirty work with MODELLER in one fell stroke.


## Tools to run MD Simulation

`pdbremix` provides a simplified cross-package interface to run a useful subset of molecular-dynamics simulations. Of course, this is in not a replacement for the full functionality of these packages.

For beginners, it is particularly useful to see how a simulation is set-up from a PDB file to a trajectory, as the shell scripts and log files of all intermediate steps are saved to file. It is easier to modify a working process than to generate one from scratch.

### Preparing Simulations from PDB

First let's grab a PDB file from the website:

	> pdbfetch 1be9

Then we can clean do some standard cleanup so that the structure exists in a unique single conformation:

	> pdbstrip 1be9.pdb

This next tool interrogates the structure for features that may affect MD simulation, highlighting steric clashes, chain-breaks (missing amino acids), disulfide bonds, incomplete and nonstandard amino acids:

	> pdbcheck 1be9.pdb

Then we generate a topology file from the PDB file:

	> pdb2sim 1be9.pdb sim AMBER11-GBSA

This will detect multiple chains, disulfide-bonds, fit hydrogen atoms to AMBER, and guess polar residue charged states. Masses, charges and bond spring parameters are generated from the AMBER99 force-field. `pdb2sim` will write a set of restart files with a common basename `sim`:

	sim.top - the toplogy file
	sim.crd - the coordinates file

The current choice of force-fields:

1. AMBER11-GBSA
2. AMBER11
3. NAMD2.8
4. GROMACS4.5

For AMBER11-GBSA, `pdb2sim` builds a topology file for implicit solvent. For the other choices, explicit solvent is used, where `pdb2sim` creates a box with 10 &Aring; padding, and fills the box with waters and counterions.

### Positional constraints

Positional constraints are very important in setting up MD simulations. `pdbremix` simplifies the application of positional restraints by using the B-factor column of PDB files to denote positional constaints, which is what NAMD does. 

To generate a PDB file for positional restraints from a set of restart files:

	> sim2pdb -b sim sim.restraint.pdb

which will generate a PDB file where all backbone atoms have been selected. You can directly edit the B-factors in the PDB file. Another option `-a` is for all protein atoms:

	> sim2pdb -a sim sim.restraint.pdb



### Running simulations

`pdbremix` provide several tools to run MD simulations where the chosen package is detected by the extensions of the restart files. 

For all packages, a robust set of simulation parameters are used, including a 1 fs time-step, and no bond-constraints on protein atoms. In explicit solvent, periodic boundary conditions are applied with PME electrostatics.

The output restart files and trajectories are written to a common basename, and an optional `-r` flag to load positional restraints:

1. Minimize your structure from `sim` restart files to `min`, using restraints defined in `sim.restraint.pdb`:

		> simmin -r sim.restraint.pdb sim min

2. MD simulation with a Langevin thermometer at 300K for 5000 fs:

		> simtemp -r restraint.pdb min temp 300 5000

3. For constant energy for 5000 fs:

		> simconst -r restraint.pdb min const 5000

This allows you to run equilibration protocols from the command-line. For instance, a prequilibration at 300K, intially a 10 ps heating of the solvent, followed by 10 ps of the system:

	> sim2pdb -b restraint.pdb sim
	> simmin -r restraint.pdb sim min
	> simtemp -r restraint.pdb min heat1 10000 300
	> simtemp heat1 heat2 10000 300
 
### Trajectory analysis

`pdbremix` provides a tool to calculate RMSD and kinetic energy for trajectories, convenience tools for viewing trajectories in viewers, and some translation tools. To use these tools, the trajectory files must have the following naming structure:

- AMBER: 
	- md.top
	- md.trj
	- md.vel.trj
- GROMACS: 
	- md.top (and md.\*itp)
	- md.gro
	- md.trr
- NAMD: 
	- md.psf
	- md.dcd
	- md.vel.dcd

These are trajectory analysis tools:

- `trajstep` displays basic parameters of a trajectory
- `trajvar` calculates energy and RMSD of trajectory

As opening trajectories in standard viewers are a pain, use these tools to open them:

- `trajvmd` display trajectory in VMD *\*recommended\**
- `trajchim` display trajectory in CHIMERA
- `trajpym` display trajectory in PYMOL *\*AMBER only\**

And some package specific tools: 

- `traj2amb` converts NAMD/GROMACS to AMBER trajectories *\*without\** solvent
- `grotrim` trim GROMACS .trr trajectory files


## Python interface to PDB structures

An important part of `pdbremix` is the design of a light API to interact with PDB structures. The data structures are designed to be easy to use with idomatic Python to do things such as select atoms. 

Other packages sometimes include a domain-specific language for atom selection, but ultimately this limits the ability for those libraries to interact with the Python ecosystem such as scipy, pandas, or numpy.


### Vector geometry library

As in any structural biology library, `pdbremix` proivdes a full-featured vector geometry library `v3`:

	from pdbremix import v3

`v3` was designed to be function-based, which allows the library to switch between a pure Python version and a numpy-dependent version. 

If you want just the python version:

	import pdbremix.v3array as v3

Or the numpy version: 

	import pdbremix.v3numpy as v3

Vectors are created and copied by the `vector` function:

	v = v3.vector() # the zero vector
	z = v3.vector(1,2,3)
	w = v3.vector(z) # a copy

Vectors are represented as arrays as they are subclassed from Python arrays or numpy arrays, and components are accessed as:

	print v[0], v[1], v[2]

All vectors functions return by value, with the one exception of `set_vector`, which changes components in place:

	v3.set_vector(v, 2, 2, 2)

Here are a set of common vector operations:

	mag(v)
	scale(v, s)
	dot(v1, v2)
	cross(v1, v2)
	norm(v)
	parallel(v, axis)
	perpendicular(v, axis)

Vectors will be used to represent coordinates/points, velocities, displacements, etc. Functions are provided to measure their geometric properties:

	distance(p1, p2)
	vec_angle(a, b)
	vec_dihedral(a, axis, c)
	dihedral(p1, p2, p3, p4)
	
	normalize_angle(angle)
	degrees(radians)
	radians(degrees)
	
	get_center(crds)
	get_width(crds)

#### Affine Transforms

We also need a representations for affine transforms, which involve a rotation and a translation. Such a transform is represented as a `matrix` that is designed to transform a vector `v`:

	v3.transform(matrix, v)

Most of the time, you would build a transform from these basic generating functions:

	v3.identity()
	v3.rotation(axis, theta)
	v3.rotation_at_center(axis, theta, center)
	v3.translation(t)
	v3.left_inverse(m)

And combine them in the correct sequence:

	c = v3.combine(a, b)

If you need to access the transform directly, transforms are represented as a 4x3 matrix, which are accesed by:

	matrix_elem(matrix, i, j, val=None)

Like any affine transform matrix, it consists of:

1. 3x3 rotational component:

		matrix_elem(m, i, j) for i=0..3, j=0..3

2. 3x1 translational component:

		matrix_elem(m, 3, i) for i=0..3 


#### Testing Vectors and Transforms

As vectors are built of floats, comparison functions are needed to test for equality within limits:
 
	is_similar_mag(a, b, small=0.0001)
	is_similar_matrix(a, b, small=0.0001)
	is_similar_vector(a, b, small=0.0001)

These functions provide random generators, useful for testing:

	random_mag() # random positive float from [0, 90]()
	random_real() # random float from [-90, 90]() 
	random_vector()
	random_rotation()
	random_matrix()


### Reading a PDB into a Soup

Here, we look at how to manipulate PDB structure. First, let's grab a PDB structure from the website using the `fetch` module:

	from pdbremix import fetch
	fetch.get_pdbs_with_http('1be9')

The main object for manipulating PDB structures is the Soup object in `pdbatoms`. We can read a Soup from a PDB file:

	from pdbremix import pdbatoms
	soup = pdbatoms.Soup("1be9.pdb")

A Soup is essentially a collection of atoms, which we can grab by:

	atoms = soup.atoms()

We can transform the soup using affine transforms from `v3`:

	displacment = v3.vector(1,0,0)
	translation = v3.translation(displacement)
	soup.transform(translation)

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

`atom.pos` and `atom.vel` are vectors defined from `v3` library. 

Atom can be indivdually transformed:

	atom.transform(translation)

To search through a bunch of atoms, you iterate through it as ... a Python list:

	hydrogens = []
	for atom in soup.atoms():
	 if atom.elem == 'H':
	   hydogens.append(atom)

Or we can leverage better Python idioms such as the `re` module:
  
	import re
	backbone_atoms = filter(
	  lambda a: re.match(r'(C|O|N|H|CA)', a.type),
	  soup.atoms())

You can search with geometric parameters:

	from pdbremix import v3
	from pdbremix import pdbatoms
	
	center = pdbatoms.get_center(soup.atoms())
	furtherest_atom = None
	furtherest_dist = 0.0
	for atom in soup.atoms():
	  d = v3.distance(center, atom.pos)
	  if d > furtherest_dist:
	    furtherest_atom = atom
	    furtherest_d = d

After, doing all this analysis, you'd probably want to save your results to a PDB file, whether to start a simulation, viewing, or interface with other programs:

	soup.write_pdb('out.pdb')

A very common strategy for saving results in PDB files, is to write residue/atom values in the B-factor column. To set the B-fcators by residue, say in:

	residue_asa = [...]

You can load into a soup by:

	soup.load_residue_bfactors(residue_asa)
	soup.write_pdb('bfactor.pdb')


### Searching through residues

The Soup object is not just a list of Atoms, it also contains a list of Residues that organize the list of Atoms. 

Following the PDB, a residue represents a collection of atoms that forms a recognisable chemical group, whether a distinct molecule in the case of solvent, ions and ligands, or an actual residue that forms a polymer, as in amino acids in a protein or a nucleic acid in DNA. This is actually quite a useful organizational definition. If we were to force the grouping of atoms in terms of molecules, then you'd be treating proteins with tens of thousands of atoms to a water molecules with 3. 

To access the list of residues:

	residues = soup.residues()

Individual residues can be in a Soup by a list index:

	first_residue = soup.residue(0)

or by a *\*residue tag\**:

	key_residue = soup.residue_by_tag('A:15').name('CA')

Each residue contains a group of atoms:

	residue = soup.residue(0)
	atoms = residue.atoms()
 
Each atom in a residue should have a unique atom type:

	ca = residue.atom("CA")

You can loop through residues and atoms. To find all atoms surrounding a residue:

	for residue in soup.residues():
	  if residue.has_atom("CA"):
	    ca = residue("CA")
	    for atom in soup.atoms():
	      if v3.distance(ca.pos, atom.pos) < 4:

To help with searching, a number of useful default lists are provided in a data module:

	from pdbremix import data
	
	data.res_name_to_char - a dictionary to convert 3 letter amino acid   
	                        names to char
	data.res_char_to_name - the reverse of the above
	data.backbone_atoms - common backbone names for PDB and MD packages
	data.radii - dictionary of atomic radii for various elements
	data.solvent_res_types - some common solvent residue names in PDB   
	                         MD packages


### Handling chains in PDB files

The Soup object does not represent chains in its data structure. In the author's experience, this adds a lot of complexity that doesn't result in any great utility. Rather, the `chain_id` is stored in every atom, and should be the same for all atoms in a residue.

The chain id's can be fetched from a Soup:

	chain_ids = soup.get_chain_ids()

And a new Soup can be extracted by chain\_id:

	chain_a = soup.get_chain('A')

If you need to deal separately with chains, store it in a Python dictionary:

	system = {}
	for chain_id in soup.get_chain_ids():
	  system[chain_id] = soup.get_chain(chain_id)

### Patching PDB structures

PDB files are complicated and messy. This is because nature is complicated and messy, and PDB files have a simple structure. To prepare PDB files for MD simulations, the PDB files might have to be heavily edited to form a single unique conformations. 

Much of this can be done with some straightforward text processing that throw away experimental information. `pdbtext` provides some of text trasnformation functions:

- `strip_lines(pdb_txt, tag_func)` - is a utility function that uses anonymous function to filter lines of texts. lines that `tag_func` return as True will be skipped.
- `strip_hydrogens(pdb_txt)` - strips out any ATOM lines that contain hydrogens
- `strip_solvent` - deletes any entries that match solvents defined in `data.solvent_res_types`
- `renumber_residues(pdb_txt) ` - cleans up the numbering of residues by renumbering residues sequentially, overwriting inserts and missing numbers
- `strip_other_nmr_models(pdb_txt) ` - NMR files, and homology models often contain many alternate but similar conformations. This function takes only the first one encountered.
- `strip_alternative_atoms(pdb_txt) ` - X-ray structures sometimes include several conformations for a residue, as it flips back and forth. This function keeps only the first conformation.

All these functions are run with:

- `clean_pdb (in_pdb, out_pdb)` - runs all the above


### Structure Analysis

Now you've happily loaded a PDB structure into Soup, you might want to do some more involved geometric analysis, with these libraries:

`asa` - calculates the accessible surface-area of every atom in list of atoms, with respect to the other atoms.

	from pdbremix import asa
	asa.calculate_asa(atoms, probe, n_sphere_point=960)

which assigns the asa to to each `atom.asa`

`volume` - calculates the volume of a list of atoms, using a standard lookup table of radii in data.radii.

	from pdbremix import volume
	volume = volume.volume(atoms, 0.5)

`rmsd` - calculates the RMSD between sets of coordinates. Two algorithms are provided, the standard SVD method when numpy is available, and the qcp algorithm of Dogulas Theobald otherwise:

	from pdbremix import rmsd

If you've extracted the positions in a list, you can calculate directly using:
 
	rmsd, transform12 = rmsd.calc_rmsd_rot(crds1, crds2)

Or for PDB structures:

	rmsd, transform12 = rmsd_of_pdbs(
	    pdbs1, pdbs2, segments1=[], segments2=[], 
	    atom_types=['CA'], transform_pdb1=None)

To use, you set up segments like this:

	segments1 = [('A:1', 'A:2'), ('B:1', 'B:10')]
	semgents2 = [('A:5', 'A:7'), ('B:2', 'B:11')]
	rmsd, transform12 = rmsd_of_pdbs(pdbs1, pdbs2, segements1, segments2)


### Making PNG of Proteins

Another common problem is generating a lot of images of proteins. `pdbremix` provides a conveninent interface to use PYMOL as a image generator library:

	from pdbremix import pymol
	pymol.make_pdb_png(
	    png, pdbs, bgcolor="white", 
	    center_res=None, top_res=None,
	    highlight_res=None, is_sticks=True,
	    is_putty=False, 
	    width=480, height=480)

This makes a `png` from a bunch of `pdbs`. The frame of viewing is defined by `center_res` and `top_res`. The default coloring is by chain, with ribbons, and sidechains as ribbons. However, if `is_putty` is false, it shows a B-factor colored worm view and the png size is determined by `width` and `height`.


## Python Interface to Molecular Dynamics

`pdbremix` provides an extensive API to run molecular-dynamics simulations from Python. `pdbremix` abstracts the particular details of running MD under different packages. 

One key abstraction in `pdbremix` is that every simulation is started with a set of restart files that share a common basename, e.g. for AMBER:

	- sim.top: topology file
	- sim.crd: coordinates/velocities file(s)

At any later point, the restart files can be obtained by the basename:

	top, crds, vels = simulate.get_restart_files('sim')

The `vels` parameter is needed for some MD packages, and is a dummy variable for other packages.

This `simulate.py` module provides a suite of functions to build restart files and run minimizations, fixed-temperature dynamics and constant-energy dynamics, and simplified the process of setting up positional restraints.

A key parameter used in these functions is `force_field` which describes the external MD-package for `pdbremix` to use:

- AMBER11-GBSA: AMBER 99 force-field with X igb surface area
- ABMER11: using the AMBER 99 force-field and waters?
- NAMD2.8: using the AMBER 99 force-field and waters
- GROMACS4.5: using the CHARM22 force-field and TIP3P waters

`pdbremix`provides a restricted set of MD simulations with proteins that is reasonably consistent across packages. They can be essentially classed into two types:

1. explicit water simulations:
	- periodic box with 10 Å padding from protein
	- waters generated to fill in the periodic box
	- Langevin thermometer
	- Nose-Hoover barometer set to 1 Atm
	- Particle-Ewald-Mesh Electrostatics
	- no bond constraints on protein
	- 1 fs timestep
2. implicit-solvent simulations:
	- Generalized Born electrostatics
	- Surface Area hydrophobic term
	- Langevin thermemeter
	- no bond constraints on protein
	- 1 fs timestep

The `force_field` of AMBER11-GBSA is implicit-solvent whilst all the others are explicit solvent.


### Restart files: topologies and coordinates

But before anything happens, you need to turn a PDB file into a set of restart files. This is carried out with the `simulate.pdb_to_top_and_crds ` function, which will take a PDB file and given the `force_field`, using the tools in each MD package:

- setup disulfide bonds
- detect charged/polar residue state
- build all hydrogens
- setup the terminii of protein chains

If the`force_field` is for explicit solvent, then the function will also:

- set up a box with 10 Angstrom padding
- fill the box with waters
- add counterions to neutralize the system

Here's an example:

	top, crds = simulate.pdb_to_top_and_crds('AMBER11', '1be9.pdb', 'sim')

This will produce the restart files for AMBER:

	- top: sim.top
	- crds: sim.crd

For GROMACS4.5, this would be:

	- top: sim.top
	- crds: sim.gro

And for NAMD2.8:

	- top: sim.psf
	- crds: sim.coor


### Positional Restraints

For a lot of MD simulation protocols, positional restraints are required. In NAMD, these are really easy to implement, simply take a PDB of the simulation, and set the Bfactors to 1 of the atoms to restrain. It's messier for AMBER and GROMACS.

`pdbremix` has abstracted positional restraints for all  3 packages to follow the NAMD methodology. You can use create such a PDB file like this:

	top, crds, vels = simulate.get_restart_files('sim.pdb')
	soup = simulate.soup_from_restart_files(top, crds, vels)
	for atom in soup.atoms():
	  if atom.type = 'CA':
	    atom.bfactor = 1.0
	  else:
	    atom.bfactor = 0.0
	soup.write_pdb('sim.rsetraint.pdb')

This `sim.restraint.pdb` can now be used for the parameter `restraint_pdb` in the following simulation functions.

### Overview of MD strategies

A typical molecular-dynamics typically involves two phases:

1. Equilibration: a prepartion stage 
2. Production: used for analysis

The equilibration can vary quite a bit for different problems. It typically consists of a mix of minimizations and constant-temperature simulations with different times and positional restraints. `pdbremix` provides three basic simulation functions:

	1. simulate.minimize
	2. simulate.langevin_thermometer
	3. simulate.constant_energy

which can be used the set up for different sequences for  equilibration, and then, for the production run.

### Minimization

One should only start a dynamics simulation with a well-minimized structure. For various reasons, the intial conformation of the PDB may place parts of the protein in a tightly compressed conformation that will burst open once a dynamics simulation starts. A minimization can find the nearest conformation where all the atoms are relaxed, and not liable to blow up. This is a better starting point for dynamics simulations.

In `pdbremix` minimizations are performed with:

	> simulate.minimize('AMBER11', 'sim', 'min')

with optional positional restraints.

### Langevin Thermometer

The most common way to run MD simulations is to hold the system artifically at a fixed temperature. `pdbremix` provides a Langevin thermometer that maintains the average velocities of the system to the target temperature by applying a stochastic force. A Gamma value of 5 is used for the stochastic force. Langevin thermometers are better than Anderson thermometers as they avoid getting trapped in certain local minima for the cost of stochasticity. This simulation is run:

	> simulate.langevin_thermometer('AMBER11', 'sim', 'min', 20000, 300)

### Constant Energy

Nevertheless, sometimes you want to see the system without the imposition of a thermometer. If the system has been pre-equilibrated to a given temperature, then you can run it a constant energy for a while. 

	> simulate.constant_energy('AMBER11', 'sim', 'min')

Unfortunately, due to the nature of numerical simulations, the integrity of the system will degrade over the length of the simulation, and the energy of the system will fluctuate.

### Examples of Equilibration Strategies

Here's an example: I wanted to do some low temperature RIP simulations. The structures had to be equilibrated to 10K, and the simulations had to run in constant energy.

1. Heating to 10K using Langevin for 10ps
2. Relaxation at constant energy for 10ps
3. Heating the relaxed conformation back to 10px

I could build a function that does all this in `pdbremix` and have restart files in a directory `equil_md`.

	simulate.minimize('sim', 'min') 
	simulate.langevin_thermometer('min', 'heat1', 1000, 10)
	simulate.constant_energy('heat1', 'const2', 1000)
	simulate.langevin_thermometer('const1', 'heat3', 1000)


### PUFF approach to steered molecular dynamics

The `pdbremix` libraries a particular powerful method of carrying out steered molecular dynamics simulations. This is the PUFF method (ref).

The idea is very simple. Steered molecular dynamics works by applying artificial forces to a system. This is normally implemented within a MD package, and requires quite a detailed setup.

However, If you run a simulation for a very short period of time, say 100 fs. Then by directly changing the velocities of the restart files, you are effectively applying forces to the system.

Since `pdbremix` possesses routines to read/write restart files into Soup objects, this can be easily done. The `simulate` module provides a function `pulse` which implements this, and the function takes as a paramter, the `pulse_fn` that performs the velocity changes to a soup.

The `force` module provides various functions that builds `pulse_fn` to carry out different types of forces:

1. Applies a Random Gas Force to the ith residue of a soup:

		make_atd_fn(i_residue, heating_temperature, backbone_atoms)
  
2. Applies a pushing force that is set to target velocity between two domains

		 make_puff_fn(
		    domain1, domain2, target_val, dt=0.1, temperature=None, 
		    is_backbone_only=False, is_first_domain_only=False, 
		    force_fname='md.puff.out')

3. Applies a pushing force that is a set acceleration:

		make_puff_acc_fn(
		    domain1, domain2, target_val, dt=0.1, temperature=None, 
		    is_backbone_only=False, is_first_domain_only=False, 
		    force_fname='md.puff.out')

4. Applies a random rotational velocity to the chi angles of the residue:

		make_rip_fn(i_res, heating_temperature)

To use here's an example:

	 top, crds, vels = simulate.get_restart_files(md)
	 n = len(simulate.soup_from_restart_files(top, crds, vels).residues())
	 pulse_fn = force.make_puff_fn([0, 1, 2], [n-3, n-2, n-1], 10.0, 0.1, 300)
	 simulate.pulse(ff, md, 'md', 2000, pulse_fn, 100)


### Reading Trajectories

Reading a trajectory in `pdbremix` is simply plugging in a Soup read from topology files to a frame reader. This is achieved from a Trajectory class defined for each package with a common interface. This interace is:

	class TrajectoryReader
	  Attributes:
	    soup
	    i_frame
	    self.n_frame
	  Method:
	    load_frame(i_frame)

You an access it through the factory function:

	import simulate.trajectory
	trj = simulate.trajectory.open_trajectory('md')

which will look for the topology and trajectories with the same base name.

To run through the frames, simply:

	trj.load_frame(52)

and then access the soup:

	atom = trj.soup.residue(0).atom("CA")
	trj.soup.write_pdb('frame.pdb')

As the `pos` and `vel` vectors of each atom are updated in place, you can safely reference atoms and residues in other variables, and their values will update as you `load_frame`:

	hydrogens = filter(lambda a: a.elem=="H", trj.soup.atoms())
	for i in range(len(trj.n_frame)):
	  trj.load_frame(-1)
	  print [h.pos for h in hydrogens]

To save a frame of a trajectory, simply copy the soup:

	trj.load_frame(0)
	first_frame_soup  = trj.soup.copy()

For convenience, there is a TrajectoryAnalysis class that simplifies the processing of trajectory analysis. To use, simply subclass TrajectoryAnalysis and set the name of the variable:

	class MyTraj(trajectory.TrajectoryAnalyzer):
	  var_name = "velocity"

Then override:

	def calulate_results(self):
	  pos = pdbatoms.get_center([a for a in self.trj.soup.atoms()])
	  return [pos]

The return list should be a list of values that will be written to the output file. Then run:

	analyze_trajectory('md', analyzer_classes=[MyTraj])

which will produce the file `md.velocity.per_frame` and `md.velocity.per_ps`, which are simply text files where each line is a list of values as returned by `calculate_results`.

[1]:	https://github.com/boscoh/pdbremix/archive/master.zip