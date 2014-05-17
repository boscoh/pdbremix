title: pdbremix documentation
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

&nbsp; &nbsp; [Zip-Package][1]

And install:

	> python setup.py install

From here, you can access unit tests and example files.

There are many wonderful tools in structural biology that have less-than-stellar interfaces. `pdbremix` wraps these tools with extra functionality. To check which tools can be access from the path:

	> checkpdbremix

If you know where the binaries are, or if you want to add exotic flags to the binaries, edit the confguration file with the `-o` flag:

	> vi `checkpdbremix -o`


## Tools to analyze PDB structures

`pdbremix` provides a bunch of tools to investigate PDB structures.

### Tools in Pure Python

Some of them can be used straight out of the box:

- `pdbfetch` fetches PDB files from the RCSB website
- `pdbheader` displays summary of PDB files
- `pdbseq` displays sequences in a PDB
- `pdbchain` extracts chains from a PDB
- `pdbcheck` checks for common defects in a PDB
- `pdbstrip` cleans up PDB for MD simulations

These tools implement standard structural biology algorithms:

- `pdbvol` calculates volume of a PDB
- `pdbasa` calculates accessible surface-area of a PDB
- `pdbrmsd` calculates RMSD between PDB files

For these algorithmic tools, you get a large speed gain if you run them  through `pypy`:

	> pdbfetch 1be9
	> pdbstrip 1be9.pdb
	> pypy `which pdbvol` 1be9.pdb

Help for the tools is available with the `-h` option. 


### Wrappers around External Tools

These following tools wrap external tools to solve some very common (and painful) use-cases in PDB analysis.

- `pdbshow` displays PDB structures in PYMOL with extras.

	PYMOL is a powerful viewer, but it's defaults leave a little to be desired. `pdbshow` runs PYMOL with some useful added functionality:

	  1. By default: colored chains, ribbons, and sidechains as sticks. 
	  2. Choice of initial viewing frame, defined by a center-residue and a top-residue. The PDB will be rotated such that the center-residue is above the center-of-mass in the middle of the screen. The top-residue will be above the center-residue.
	  3. Color by B-factor using a red-white scale, with limits defined by options.
	  4. Solvent molecules can be removed, which is specifically for MD frames that contain too many waters.

- `pdboverlay` display homologous PDB files using MAFFT, THESEUS and PYMOL.

	One of the most beautiful results of structural biology is the structural alignment of homologous proteins. `pdboverlay` performs this complex process in one easy step starting from PDB structures:

	1. Write fasta sequences from PDB.
	2. Align sequences with MAFTT to find homologous regions.
	3. Structurally align homologous regions with THESEUS.
	4. Display structurally-aligned PDBs using special PYMOL script.

- `pdbinsert` fill gaps in PDB with MODELLER

	Gaps in PDB structures cause terrible problems in MD simulations. The standard tool to patch gaps is MODELLER, which requires a ton of boilerplate. `pdbinsert` does all the dirty work with MODELLER in one fell stroke.


## Tools to run MD Simulation

`pdbremix` provides a simplified interface to run molecular-dynamics. Whilst this is not a replacement for understanding each of the MD packages, the scenarios provided here cover a reasonable range of cases.

For beginners, it is useful to see how a basic simulation is set-up from a PDB file to a trajectory, as all intermediate files and shell scripts of the commands are are saved to file as well as full logging for all packages. It is much easier to modify a complete and working process than to start from scratch.

### Preparing Simulations from PDB

First let's grab a PDB file from the website::

	> pdbfetch 1be9

Then we can clean it up into a single conformation ready for MD:

	> pdbstrip 1be9.pdb

Hopefully everything is good:

	> pdbcheck 1be9.pdb

Then we generate a topology file from the PDB file:

	> pdb2sim 1be9.pdb sim AMBER11-GBSA

This will detect multiple chains, disulfide-bonds, fit hydrogen atoms to AMBER, and guess polar residue charged states. Masses, charges and bond spring parameters are generated from the AMBER99 force-field. The restart is a set of restart files with a common basename:

1. sim.top - the toplogy file
2. sim.crd - the coordinates

Now there are several choice of packages/force-fields:

1. AMBER11-GBSA
2. AMBER11
3. NAMD2.8
4. GROMACS4.5

The first one builds a topology file for implicit solvent. The rest builds explict solvent. For explicit solvent, `pdb2sim` also creates a box with 10 &Angs; padding, and fills the box with waters and coutnerions.

### Positional constraints

Positional constraints are very important in setting up MD simulations. `pdbremix` simplifies the application of positional restraints by using the B-factor column of PDB files to denote positional constaints, which is what NAMD does. If you have a set of restart files, then run:

	> sim2pdb -b sim sim.restraint.pdb

which will generate a PDB file where all backbone atoms have been set to B-factor=-1. This will be used to apply positional restraints. Another option is:

	> sim2pdb -a sim sim.restraint.pdb

Of course you can edit your own B-factors in the PDB file.

### Running simulations

Typically you want to run a simulation at a given temperature, so you want to start the simulation that has been heated carefully to that temperature without suffering any kind of heating artefact. 

`pdbremix` several wrappers to run MD simulations. The MD package is detected by the extensions of the restart files. A robust set of simulation parameters are chosen for you: no bond constraints on protein and a 1 fs time-step. In explicit solvent: also periodic boundary conditions and PME electrostatics. The simulation wrappers are run (with the optional `-r` restraint flag:

1. Tinimize your structure from `sim` restart files to `min`, using restraints defined in `sim.restraint.pdb`:

		> simmin -r sim.restraint.pdb sim min

2. MD simulation with a Langevin thermometer at 300K for 5000 fs:

		> simtemp -r restraint.pdb min temp 300 5000

3. For constant energy for 5000 fs:

		> simconst -r restraint.pdb min const 5000

You can then study the trajectories with the tools in the next section.


### Trajectory analysis

Trajectory analysis is probably best done with dedicated scripts, but here's some scripts for some quick analysis. I found particular useful ones are scripts that opens a trajectory in a viewer. 

Before you can use this scripts you have to make sure the trajectory files share a common basename:

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

Once named properly, these tools can be used:

- `trajstep` extracts basic parameters of a trajectory
- `trajvar` calculates energy and RMSD of trajectory

Use these tools to display trajectories: 

- `trajvmd` display trajectory in VMD *\*recommended\**
- `trajchim` display trajectory in CHIMERA
- `trajpym` display trajectory in PYMOL *\*AMBER only\**

And some package specific tools: 

- `traj2amb` converts NAMD/GROMACS to AMBER trajectories **\*without\*** solvent
- `grotrim` trim GROMACS .trr trajectory files


## Python interface to PDB structures

An important part of `pdbremix` is the design of a light API to interact with PDB structures. The data structures are meant to be as small as possible that can provide most reasonable functionality, hence the inclusion of a fully-featured vector library.

Some packages provide atom selecton using a domain-specific language. Whilst this may simplify things at the beginning, this makes it difficult to interact with the full power of Python. 

Here we prefer to use as much idomatic Python to select atoms, and to analyze geometry. By staying in Python, this makes it easier to write code to interact with other Python libraries such as sci.py or pandas, or numpy.


### Vector geometry library

As in any structural biology library, we provide a vector geometry library, which is called `v3`:

	from pdbremix import v3

`v3` was designed to be function-based, with vector and transform objects subclassed from arrays. This allows the library to easily switch between a pure Python version and a numpy-dependent version. If you want just the python version:

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

Here are a set of common vector operationsnthat returns by value:

	mag(v)
	scale(v, s)
	dot(v1, v2)
	cross(v1, v2)
	norm(v)
	parallel(v, axis)
	perpendicular(v, axis)

Vectors will be used to represent coordinates/points, velocities, displacements etc. Functions are provided to measure their geometric properties:

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

Finally, we introduce functions to test similarity for vectors and transforms as floats are not exact and require tolerance for comparison:
 
	is_similar_mag(a, b, small=0.0001)
	is_similar_matrix(a, b, small=0.0001)
	is_similar_vector(a, b, small=0.0001)

And a bunch of random geometric object generators, used for testing:  

	random_mag() # random positive float from [0, 90]()
	random_real() # random float from [-90, 90]() 
	random_vector()
	random_rotation()
	random_matrix()


### Reading a PDB into a Soup

Here, we look at how to manipulate PDB structure. First, let's grab a PDB structure from the website using the `fetch` module:

	from pdbremix import fetch
	fetch.get_pdbs_with_http(['1be9'])

The main object for manipulating PDB structures is the Soup object in `pdbatoms`. We can read a Soup from a PDB file:

	from pdbremix import pdbatoms
	soup = pdbatoms.Soup("1be9.pdb")

A Soup is essentially a collection of atoms, which we can grab by:

	atoms = soup.atoms()

An atom has attributes:

  - pos (v3.vector)
  - vel (v3.vector)
  - mass (float)
  - charge (float)
  - type (str)
  - element (str)
  - num (int)
  - chain\_id (str)
  - res\_type (str)
  - res\_num (str)
  - res\_insert (str)
  - bfactor (float)
  - occupancy (float)
  - alt\_conform (str)
  - is\_hetatm (bool)

Since `atom.pos` and `atom.vel` are vectors, these can be manipulated by functions from the `v3` library. Atom contains one special method `transform` that handle affine transforms:

	displacment = v3.vector(1,0,0)
	translation = v3.translation(displacement)
	atom.transform(translation)

Given that a soup is essentially a list of atoms, you can transform the entire soup: 

	soup.transform(translation)

To search through a bunch of atoms, you iterate through it as ... a Python list:

	hydrogens = []
	for atom in soup.atoms():
	 if atom.elem == 'H':
	   hydogens.append(atom)

Vectors play a natural part, you just have to remember to call the `v3` library and refer to `pos` and `vel`

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

- Use the Python `re` module to string analysis


### Searching through residues

The Soup object is not simply a list of Atoms, it also contains a list of Residues that organize the list of Atoms. 

A residue in this case represents a collection of atoms that forms a recognisable chemical group, whether a distinct molecule in the case of solvent, ions and ligands, or an actual residue that forms a polymer, as in amino acids in a protein or a nucleic acid in DNA. 

This definition follow the PDB is actually quite a useful organizational definition. If we were to force the grouping of atoms in terms of molecules, then you'd be treating proteins with tens of thousands of atoms to water molecules with 3. 

To access the list of residues:

	residues = soup.residues()

One heuristic of residues is that each atom in a residue has a unique atom type. Then one can look up atoms by residue and atom type:

- residue atoms

	 ca = residue.atom('CA')

Individual residues can be in a Soup by a list index:

	first_residue = soup.residue(0)

But a more useful one is by a *residue tag*:

	key_residue = soup.residue_by_tag('A:15').name('CA')

You can loop through residues and atoms. To find all atoms surrounding a residue:

	for residue in soup.residues():
	  if residue.has_atom("CA"):
	    ca = residue("CA")
	    for atom in soup.atoms():
	      if v3.distance(ca.pos, atom.pos) < 4:

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

PDB files are complicated and messy. This is because nature is complicated and messy, and PDB files are relatively simple. In order to perform analysis on them such as with `pdbatoms` and prepare them for MD simulations.

There are some simple text processing steps that throw away a lot of the experimental information so that we can end up with one single, largely intact conformation to analyze. `pdbtext` provides some of these functions:

- `strip_lines` - is a utility function that uses anonymous function to filter lines of texts
- `strip_hydrogens` - strips out any ATOM lines that contain hydrogens
- `strip_solvent` - deletes any entries that match solvents defined in `data.solvent_res_types`
- `renumber_residues` - cleans up the numbering of residues by renumbering residues sequentially, overwriting inserts and missing numbers
- `strip_other_nmr_models` - NMR files, and homology models often contain many alternate but similar conformations. This function takes only the first one encountered.
- `strip_alternative_atoms` - X-ray structures sometimes include several conformations for a residue, as it flips back and forth. This function keeps only the first conformation.
- `clean_pdb` - runs all the above 


### Structure Analysis

Now you've happily loaded a PDB structure into Soup, you might want to do some analysis.

We've provided libraries to carry out some common analysis.

`asa` - calculates the accessible surface-area of every atom in list of atoms, with respect to the other atoms.

`volume` - calculates the volume of a list of atoms, using a standard lookup table of radii in data.radii.

`rmsd` - calculates the RMSD between sets of coordinates. These can be extracted from soups, or PDB files. Two algorithms are provided, the qcp algorithm of Dogulas Theobald, or the stand SVD method using numpy.

`protein` - contains a few common protein analysis

- moving things around

Finally, you want to save your results. This is simply:

	soup.write_pdb('out.pdb')

However, you might also want to save any atom/residue analysis. Most people do this  by saving the values in the B-factor of a column. If you have a list of values, such as:

	residue_asa = [...]

You can load into a soup by:

	soup.load_residue_bfactors(residue_asa)

And if you have atoms, then do a list comprehension:

	atom_asa_list = [...]
	[atom.bfactor = asa 
	 for atom, asa in zip(soup.atoms, atom_asa_list)]


### Making PNG of Proteins

Apart from making PDB structures, another common and painful problem is to generate images of proteins.

`pymol` provides an interface to pymol to generate images of proteins in an automated way.

The very first thing that one must do is to orientate the protein to a specific frame of reference. This is very hard to do in PYMOL. 

`pdbremix` provides a simple way to get at most frames of reference. A frame of reference is simply defined by choosing a center-res and a top-res. 

The routines in `pymol` orientates the structure so that the center-res is place in the middle of the screen, directly above the center-of-mass, and rotates around this axis such that the top-res is found above the center-res.

	make_pdb_png(
	    png, pdbs, bgcolor="white", 
	     center_res=None, top_res=None,
	    highlight_res=None, is_sticks=True,
	     is_putty=False, 
	     width=480, height=480)

This shows a protein in ribbon conformation, with a highlight\_res.



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

But before anything happens, you need to turn a PDB file into a set of restart files. This is carried out with the `simulate.pdb_to_topology` function, which will take a PDB file and given the `force_field`, using the tools in each MD package:

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
	2. langevin\_thermometer
	3. constant\_energy

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

The library was originally written to do PUFF steered-molecular dynamics simulation that uses a pulsed force application. This method does not require the underlying MD package to support the method and can be carried out by manipulating restart files. The utilities to do this are:


### Reading Trajectories

Reading a trajectory in `pdbremix` is simply plugging in a Soup read from topology files to a frame reader. This is achieved from a Trajectory class defined for each package with a common interface. This interace is

	class TrajectoryReader
	  Attributes
	    soup
	    i\_frame
	    self.n\_frame
	  Method
	    __init__
	    load\_frame(i\_frame)

You an access it through the factory function:

	import simulate.trajectory
	trj = simulate.trajectory.open_trajectory('md')

which will look for the topology and trajectories with the sam base name.

To run through the frames, simply:

	trj.load_frame(52)

and then access the soup:

	atom = trj.soup.residue(10).atom("CA")

This soup is a sub-class used above for PDB analysis, and you do the same stuff here:

	trj.soup.write_pdb('frame.pdb')

As the `pos` and `v`el\` vectors of each atom is updated, you can reference them add access them individually.
\`
	hydrogens = [a for a in trj.soup.atoms() if  a.elem="H"]
	 trj.load_frame(-1)
	print [h.pos for h in hydrogens]

or make a copy of the soup to do destructive analysis:

	soup_copy  = trj.soup.copy()

Now you use the entire apparatus of the vector geometry and PDB analysis on this copy of the frame.

In `simulate.trajectory` there is a TrajectoryAnalysis class that simplifies the processing of structural parameters of a trajectory. To use, simply subclass TrajectoryAnalysis and override the function

	calulate_results()

Set the name of the variable in your file by adding a class name var\_name to the class.

that will return a list of results per frame of the trajectory. The function has access to `self.trj.soup` to process whatever parameters you want. As well, since this an object, persistent data can be preserved between frames by writing to `self`. The values you returned will be saved to a filename with the var\_name defined in the calss.

Then run analyze\_trajectory and pass in your TrajectoryAnalysis as a parameter to run the analysis.


[1]:	https://github.com/boscoh/pdbremix/archive/master.zip