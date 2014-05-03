# PDBREMIX

# Introduction

PDBSTRUCT2 is a python library for computational structural biology. It provides a unified interface to manipulate PDB structures, setup basic molecular-dynamics simulations, and analyse the resultant trajectories. 

# Standalone Command-line Utilities

For analysis and manipulation of PDB structures:

- `pdbheader` displays header information in table format
- `pdbvol` calculates the volume of the molecule
- `pdbasa` calculates the solvent-accessible surface-area
- `pdbfetch` fetches PDB files from the RCSB website
- `pdbrmsd` calculates the RMSD between PDB structures
- `pdbseq` displays the sequence actually in the PDB
- `pdbchain` extracts individual chains from a PDB file
- `pdbcheck` checks for defects that will prevent MD simulation

# Wrappers around external tools

- `pdbshow` displays PDB files with PYMOL, with extra options such as residue centring and b-factor colouring
- `pdboverlay` uses MAFFT and THESEUS to align homologous proteins and displays them with PYMOL
- `pdbinsert` uses MODELLER to build loops for gaps in a PDB structure

# Vector3d

- numerics
- function based, not object based
- allows wrapping around other libraries
- not operator overloading, with matrices, it's too easy to get confused
- allows non-numpy library, so pypy
  pypy `which pdbasa` 1be9.pdb

# reading pdb structures

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
http://web.mit.edu/vmd_v1.9.1/namd-tutorial-unix.pdf
http://bionano.physics.illinois.edu/Tutorials/ssbTutorial.pdf