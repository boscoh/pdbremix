# PDBREMIX

`pdbremix` is a python library for computational structural biology. It provides a unified interface to manipulate PDB structures, setup basic molecular-dynamics simulations, and analyse the resultant trajectories. 

The library attempts to provide the lightest API for the given functionality. For instance, even `numpy` dependency is optional. Trajectory readers are written in pure Python.

## Command-line Utilities

To use, install using pip?

Then run the command with the `-h` option 

Run with `pypy` option

### Standalone tools

`pdbremix` provides a bunch of structural biology algorithms, all written in standard Python:

- `pdbfetch` fetches PDB files from the RCSB website
- `pdbheader` displays header information in table format
- `pdbvol` calculates the volume
- `pdbasa` calculates the solvent-accessible surface-area
- `pdbseq` displays the sequence actually in the PDB
- `pdbrmsd` calculates the RMSD between two PDB structures
- `pdbchain` extracts individual chains from a PDB file
- `pdbcheck` checks common defects that affect MD simulations

### Wrappers around external tools

There are many wonderful tools for computational structural biology that have less-than-stellar command-line interfaces. `pdbremix` can be used to wrap these tools with a friendlier interface, and add some useful extra functionality.

First though, to use these tools, you have to tell `pdbremix` where the binaries are. **HOW SHOULD THIS BE DONE?**

These command-line tools provide a much better interface to commonly used tools such as PYMOL, MODELLER and THESEUS:

- `pdbshow` wrapper for PYMOL to display PDB with useful defaults and extras such as residue centring and b-factor colouring
- `pdboverlay` wraps MAFFT, THESEUS and PYMOL to display homologous proteins
- `pdbinsert` wraps MODELLER to build loops for gaps in a PDB structure

### Trajectory analysis tools

`pdbremix` provides a simplified interface to run molecular-dynamics on 3 standard MD packages: AMBER11+, GROMACS4.5+ and NAMD2.8+.

In order to this, `pdbremix` makes a strong assumption that all associated files for an MD trajectory *has a common base name*. Typically an MD simulation will include:

1. topology file
2. initial coordinate/velocity file(s)
3. trajectory file(s)
4. restart coordinate/velocity file(s)

These files must have a common base name, say `md`.If the filenames do not correspond to this naming pattern, the `pdbremix` tools will not recognise the simulation and trajectory. 

This is the path to sanity. 

In AMBER, an example set of files would be:

1. md.top
2. md.in.crd or md.in.rst
3. md.trj and md.vel.trj
4. md.rst

In GROMACS:

1. md.top and associated md.\*.itp files
2. md.in.gro
3. md.trr and md.vel.trr
4. md.gro

In NAMD:

1. md.psf
2. md.in.coor and md.in.vel
3. md.dcd and md.vel.dcd
4. md.coor and md.vel

Once named properly, these `pdbremix` tools can be used to analyse trajectories:

- `trajvar` calculates kinetic energy and RMSD on trajectories
- `md2pdb` converts MD restart files into PDB
- `trajstep` gets basic MD parameters from a trajectory
- `traj2amb` converts NAMD and GROMACS trajectories to AMBER trajectories _without_ solvent

- `trajvmd` wrapper to open trajectories in VMD *recommended*
- `trajchim` wrapper to open trajectories in CHIMERA
- `trajpym` wrapper to open trajectories in PYMOL, only works with AMBER trajectories
- `grotrim` wrapper to trim GROMACS .trr trajectory files


### PUFF steered molecular dynamics

The library was originally written to do PUFF steered-molecular dynamics simulation that uses a pulsed force application. This method does not require the underlying MD package to support the method and can be carried out by manipulating restart files. The utilities to do this are:

- `puff` runs a PUFF steered molecular-dynamics from a special configuration file
- `puffshow` displays residues that are pulled in the PUFF configuration files

## Modules in pdbremix

asa.py
data.py
fetch.py
force.py
gromacs.py
namd.py
pdbatoms.py
pdbtext.py
protein.py
pymol.py
rmsd.py
simulate.py
traj.py
util.py
v3.py
v3array.py
v3numpy.py
volume.py

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