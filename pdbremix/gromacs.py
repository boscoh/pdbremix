#!/usr/bin/env python
# encoding: utf-8
"""
gromacs.py

Created by Bosco Ho on 2009-10-13.
Copyright (c) 2009 __MyCompanyName__. All rights reserved.
"""

import os
import shutil
import glob
import xdrlib

import util
import pdbatoms
import v3
import data


# mdrun = "mpiexec -np 8 /home/bosco/Packages/gromacs-4.0.7/bin/mdrun"


# GROMACS files
# ---------------------------------------------------
# compatible force-fields:
#  - GROMACS4.0.7
#  - GROMACS4.5.5
#  - GROMACS4.5.4
#
# We will always solvate, use periodic 
# boundaries and Particle Ewald Mesh.
3
# topology: sim.top
# coordinates/velocities/restart : sim.gro
#
# positions: nanometers
# velocities: nanometers/picosecond


# Routines to read topology and restart files

def expand_restart_files(name):
  top = os.path.abspath(name + '.top')
  crds = os.path.abspath(name + '.gro')
  vels = ''
  return top, crds, vels


def get_restart_files(name):
  top = os.path.abspath(name + '.top')
  crds = os.path.abspath(name + '.gro')
  util.check_files(top, crds)
  vels = ''
  return top, crds, vels


def read_top(top, chain=""):
  util.check_output(top)
  lines = open(top).readlines()

  atoms = []
  is_chain_topologies = False
  top_dir = os.path.dirname(top)
  for l in lines:
    if not is_chain_topologies:
      if 'chain topologies' in l:
        is_chain_topologies = True
      continue
    if l.startswith("#include"):
      itp = l.split()[1][1:-1]
      itp = os.path.join(top_dir, itp)
      if os.path.isfile(itp):
        full_chain_name = os.path.splitext(itp)[0]
        chain = full_chain_name.split('_')[-1]
        these_atoms = read_top(itp, chain)
        atoms.extend(these_atoms)
    if l.startswith(";"):
      break

  is_atoms = False
  qtot = None
  for l in lines:
    if not is_atoms:
      if '[ atoms ]' in l:
        is_atoms = True
      continue
    if l.startswith('['):
      break
    if l.startswith(";"):
      continue    
    if not l.strip():
      continue
    words = l.split()
    n = int(words[0])
    res_num = int(words[2])
    res_type = words[3]
    q = float(words[6])
    mass = float(words[7])
    atoms.append((mass, q, chain))

  return atoms
  

def read_gro(gro):
  util.check_output(gro)
  positions, velocities, names, nums = [], [], [], []
  lines = open(gro).readlines()
  for l in lines[2:-1]:
    words = l.split()
    if len(l) <  43:
      continue
    positions.append([10.*float(w) for w in l[20:44].split()])
    if len(l) < 63:
      continue
    velocities.append([10.*float(w) for w in l[44:68].split()])
    names.append(l[9:15].strip())
    nums.append(int(l[15:20]))
  box = [float(w) for w in lines[-1].split()]
  return positions, velocities, box, names


def AtomFromGroLine(line):
  """
  Returns an Atom object from an atom line in a pdb file.
  Converts from nm/ps to angstroms/ps
  """
  atom = pdbatoms.Atom()
  atom.res_num = int(line[0:5])
  atom.res_type = line[5:8]
  atom.type = line[10:15].strip(" ")
  atom.element = pdbatoms.guess_element(
      atom.res_type, line[12:15])
  atom.num = int(line[15:20])
  x = 10.0*float(line[20:28])
  y = 10.0*float(line[28:36])
  z = 10.0*float(line[36:44])
  v3.set_vector(atom.pos, x, y, z)
  if len(line) > 62:
    x = 10.0*float(line[44:52])
    y = 10.0*float(line[52:60])
    z = 10.0*float(line[60:68])
    v3.set_vector(atom.vel, x, y, z)
  return atom


def convert_to_pdb_atom_names(soup):
  for res in soup.residues():
    if res.type == "ILE":
      if res.has_atom("CD"):
        res.atom("CD").type = "CD1"
    if res.has_atom("OC2"):
      res.atom("OC2").type = "OXT"
    if res.has_atom("OC1"):
      res.atom("OC1").type = "O"


def convert_to_gromacs_atom_names(soup):
  for res in soup.residues():
    if res.type == "ILE":
      if res.has_atom("CD1"):
        res.atom("CD1").type = "CD"
    if res.has_atom("OXT"):
      res.atom("OXT").type = "OC2"
      if res.has_atom("O"):
        res.atom("O").type = "OC1"


def soup_from_top_gro(top, gro, skip_solvent=True):
  """
  Returns a Polymer object built from GROMACS restart files.
  Will read all atoms until the solvent. 
  """
  util.check_output(top)
  util.check_output(gro)
  atoms = []
  temp_list = []
  lines = open(gro, 'r').readlines()
  remaining_text = ""
  n_remaining_text = 0
  for i_line, line in enumerate(lines[2:-1]):
    atom = AtomFromGroLine(line)
    if skip_solvent and atom.res_type == "SOL":
      remaining_text = "".join(lines[i_line+2:-1])
      n_remaining_text = len(lines[i_line+2:-1])
      break
    atoms.append(atom)
  box = [float(w) for w in lines[-1].split()]
  top_atoms = read_top(top)
  for a, (mass, qtot, chain) in zip(atoms, top_atoms):
    a.mass = mass # mass
    a.chain_id = chain
  soup = pdbatoms.Polymer()
  curr_res_num = -1
  for a in atoms:
    if curr_res_num != a.res_num:
      res = pdbatoms.Residue(
          a.res_type, a.chain_id, a.res_num)
      soup.append_residue(res.copy())
      curr_res_num = a.res_num
    soup.insert_atom(-1, a)
  convert_to_pdb_atom_names(soup)
  soup.box = box
  soup.remaining_text = remaining_text
  soup.n_remaining_text = n_remaining_text
  return soup


def soup_from_restart_files(top, crds, vels, skip_solvent=True):
  return soup_from_top_gro(top, crds, skip_solvent)


def convert_restart_to_pdb(md_name, pdb):
  top, crds, vels = get_restart_files(md_name)
  soup = soup_from_restart_files(top, crds, vels)
  soup.write_pdb(pdb)


def write_soup_to_gro(in_soup, gro):
  soup = in_soup.copy()
  convert_to_gromacs_atom_names(soup)
  f = open(gro, 'w')
  f.write("Generated by gromacs.py\n")
  atoms = soup.atoms()
  n_atom = len(atoms) + soup.n_remaining_text
  f.write(str(n_atom) + '\n')
  for a in atoms:
    if a.res_type == "ILE" and a.type == "CD1":
      a.type = "CD"
    # GRO doesn't care about actual numbering,
    # will wrap when the columns are fill
    res_num = a.res_num % 100000
    atom_num = a.num % 100000
    s = "%5d%-5s%5s%5d%8.3f%8.3f%8.3f%8.4f%8.4f%8.4f\n" % \
        (res_num, a.res_type, a.type, atom_num, 
         a.pos[0]*0.1, a.pos[1]*0.1, a.pos[2]*0.1,
         a.vel[0]*0.1, a.vel[1]*0.1, a.vel[2]*0.1)
    f.write(s)
  if soup.remaining_text:
    f.write(soup.remaining_text)
  f.write("%10.5f%10.5f%10.5f\n" % (soup.box[0], soup.box[1], soup.box[2]))
  f.close()
  

def write_soup_to_crds_and_vels(soup, name):
  write_soup_to_gro(soup, name + '.gro')
  return name + '.gro', ''
  

def convert_to_pdb(gro, pdb):
  data.binary('editconf', '-f %s -o %s' % (editconf, gro, pdb), gro)
  # fix the naming conflict with ILE CD
  lines = open(pdb, 'r')
  new_lines = []
  for l in lines:
    res = l[17:20]
    atom = l[13:16]
    if res == "ILE" and atom == "CD ":
      new_l = l[:13] + "CD1" + l[16:]
    else:
      new_l = l
    new_lines.append(new_l)
  open(pdb, 'w').write(''.join(new_lines))


# Routines to generate topology and gro files

ions_mdp = """
; ions.mdp - used as input into grompp to generate ions.tpr
; Parameters describing what to do, when to stop and what to save
integrator  = steep   ; Algorithm (steep = steepest descent minimization)
emtol       = 1000.0  ; Stop minimization when the maximum force < 1000.0 kJ/mol/nm
emstep      = 0.01    ; Energy step size
nsteps      = 50000   ; Maximum number of (minimization) steps to perform

; Parameters describing how to find the neighbors of each atom and how to calculate the interactions
nstlist     = 1       ; Frequency to update the neighbor list and long range forces
ns_type     = grid    ; Method to determine neighbor list (simple, grid)
rlist       = 1.0     ; Cut-off for making neighbor list (short range forces)
coulombtype = PME     ; Treatment of long range electrostatic interactions
rcoulomb    = 1.0     ; Short-range electrostatic cut-off
rvdw        = 1.0     ; Short-range Van der Waals cut-off
pbc         = xyz     ; Periodic Boundary Conditions (xyz/no)
"""


def neutralize_system_with_salt(
    in_top, in_gro, out_name, force_field):
  atoms = read_top(in_top)
  qtot = sum([q for mass, q, chain in atoms])
  counter_ion_charge = -int(qtot)
  if counter_ion_charge == 0:
    shutil.copy(in_gro, out_name + '.gro')
    return

  in_mdp = out_name + '.salt.grompp.mdp'
  open(in_mdp, 'w').write(ions_mdp)
  top = out_name + '.top'
  if in_top != top:
    shutil.copy(in_top, top)
  tpr = out_name + '.salt.tpr'
  out_mdp = out_name + '.mdp'
  data.binary(
      'grompp',
      '-f %s -po %s -c %s -p %s -o %s' \
            % (in_mdp, out_mdp, in_gro, top, tpr),
      out_name + '.salt.grompp')
  util.check_files(tpr)

  gro = out_name + '.gro'
  log = out_name + '.salt.genion.log'

  # Genion requires user input to select SOL for ion replacement
  input_fname = out_name + '.salt.genion.in'
  open(input_fname, 'w').write('SOL') # pick solvent

  charge_str = ""
  if 'GROMACS4.5' in force_field:
    charge_str = " -pname NA -nname CL "
  elif 'GROMACS4.0' in force_field:
    charge_str = " -pname NA+ -nname CL- "
  if not charge_str:
    raise ValueError, "Cannot recognize force_field " + force_field
  if counter_ion_charge > 0:
    charge_str += " -np %d " % counter_ion_charge
  else:
    charge_str += " -nn %d " % abs(counter_ion_charge)
  data.binary(
      'genion', 
      '-g %s -s %s -o %s -p %s -neutral %s' % \
          (log, tpr, gro, top, charge_str),
      out_name + '.salt.genion',
      input_fname)
  util.check_files(gro)


def pdb_to_topology(
    force_field, pdb, name, 
    solvent_buffer=10):

  util.check_files(pdb)
  full_pdb = os.path.abspath(pdb)

  solv_dir = name + '.solvate'
  save_dir = os.getcwd()
  util.goto_dir(solv_dir)

  pdb_copy = os.path.basename(pdb)
  if not os.path.isfile(pdb_copy):
    shutil.copy(full_pdb, pdb_copy)

  pdb2gmx_gro = name + '.pdb2gmx.gro'
  top = name + '.top'
  itp = name + '_posre.itp'

  if 'GROMACS4.5' in force_field:
    ff = 'amber99' # much better than the default gromos43a1
  elif 'GROMACS4.0' in force_field:
    ff = 'G43a1' 
  else:
    raise ValueError, "Couldn't work out pdb2gmx for " + force_field
     
  args = '-ignh -ff %s -water spc -missing -f %s -o %s -p %s -i %s' \
          % (ff, pdb, pdb2gmx_gro, top, itp)
  data.binary('pdb2gmx', args, name+'.pdb2gmx')
  util.check_files(pdb2gmx_gro)

  box_gro = name + '.box.gro'
  solvent_buffer_in_nm = solvent_buffer/10.0 
  data.binary(
      'editconf', 
      '-f %s -o %s -c -d %f -bt cubic' \
          % (pdb2gmx_gro, box_gro, solvent_buffer_in_nm),
      name+'.box')
  util.check_files(box_gro)

  solvated_gro = name + '.solvated.gro'
  data.binary(
      'genbox',
      '-cp %s -cs spc216.gro -o %s -p %s' \
          % (box_gro, solvated_gro, top),
       '%s.solvated' % name)
  util.check_files(solvated_gro)

  gro = name + '.gro'
  neutralize_system_with_salt(top, solvated_gro, name, force_field)
  util.check_files(gro)

  convert_restart_to_pdb(name, name+'.pdb')

  fnames = util.re_glob(
      '*', os.path.basename(name) + r'[^\.]*\.(pdb|itp|gro|mdp|top)$')
  for fname in fnames:
    shutil.copy(fname, save_dir)
  delete_backup_files(name)

  os.chdir(save_dir)

  return top, gro


# simulation routines


def delete_backup_files(tag):
  util.clean_fname(*util.re_glob('*', '^#' + tag))
    

minimization_parms = { 
  'topology' : 'in.top', 
  'input_crds' : 'in.gro', 
  'output_name' : 'min', 
  'force_field': 'GROMACS',
  'restraint_pdb': '',
  'n_step_minimization' : 100, 
} 

min_mdp = """
; template .mdp file used as input into grompp to generate energy minimization for mdrun

; Parameters describing what to do, when to stop and what to save
integrator  = steep    ; Algorithm (steep = steepest descent minimization)
emtol       = 1000.0   ; Stop minimization when the maximum force < 1000.0 kJ/mol/nm
emstep      = 0.01     ; Energy step size
nsteps      = %(n_step_minimization)s    ; Maximum number of (minimization) steps to perform

; Parameters describing how to find the neighbors of each atom and how to calculate the interactions
nstlist     = 1        ; Frequency to update the neighbor list and long range forces
ns_type     = grid     ; Method to determine neighbor list (simple, grid)
rlist       = 1.0      ; Cut-off for making neighbor list (short range forces)
coulombtype = PME      ; Treatment of long range electrostatic interactions
rcoulomb    = 1.0      ; Short-range electrostatic cut-off
rvdw        = 1.0      ; Short-range Van der Waals cut-off
pbc         = xyz      ; 3-dimensional periodic boundary conditions (xyz|no)
"""

constant_energy_parms = { 
  'topology' : 'in.top', 
  'input_crds' : 'in.gro', 
  'output_name' : 'md', 
  'force_field': 'GROMACS',
  'solvent_state': 2,
  'surface_area': 1,
  'restraint_pdb': '',
  'n_step_per_snapshot' : 50, 
  'n_step_dynamics' : 1000, 
  'temp_initial': 0.0, # ignored if it is 0.0
} 

langevin_thermometer_parms = { 
  'topology' : 'in.top', 
  'input_crds' : 'in.gro', 
  'output_name' : 'md', 
  'force_field': 'GROMACS',
  'cutoff': 16.0,  # non-bonded cutoff
  'restraint_pdb': '',
  'random_seed' : 2342, 
  'temp_thermometer' : 300.0, 
  'temp_initial': 0.0, # ignored if it is 0.0
  'n_step_per_snapshot' : 50, 
  'n_step_dynamics' : 1000, 
} 

dynamics_mdp = """
; Run parameters
integrator      = md            ; leap-frog integrator
nsteps          = %(n_step_dynamics)s  ; time = n_step_dynamics*dt
dt              = 0.001         ; time-step in fs
; Output control
nstxout         = %(n_step_per_snapshot)s  ; save coordinates every 0.05 ps
nstvout         = %(n_step_per_snapshot)s  ; save velocities every 0.05 ps
nstenergy       = %(n_step_per_snapshot)s  ; save energies every 0.05 ps
nstlog          = %(n_step_per_snapshot)s  ; update log file every 0.05 ps
; Bond parameters
continuation    = yes           ; first dynamics run
constraints     = hbonds        ; bonds from heavy atom to H, constrained
; constraints     = all-bonds     ; all bonds (even heavy atom-H bonds) constrained
; constraint_algorithm = lincs    ; holonomic constraints 
; lincs_iter      = 1             ; accuracy of LINCS
; lincs_order     = 4             ; also related to accuracy
; Neighborsearching
ns_type         = grid          ; search neighboring grid cels
nstlist         = 5             ; 10 fs
rlist           = 1.0           ; short-range neighborlist cutoff (in nm)
rcoulomb        = 1.0           ; short-range electrostatic cutoff (in nm)
rvdw            = 1.0           ; short-range van der Waals cutoff (in nm)
; Periodic boundary conditions
pbc             = xyz           ; 3-dimensional periodic boundary conditions (xyz|no)
; Electrostatics
; coulombtype     = PME           ; Particle Mesh Ewald for long-range electrostatics
; pme_order       = 4             ; cubic interpolation
; fourierspacing  = 0.16          ; grid spacing for FFT
; Pressure coupling is on
pcoupl          = Parrinello-Rahman   ; Pressure coupling on in NPT
pcoupltype      = isotropic     ; uniform scaling of box vectors
tau_p           = 2.0           ; time constant, in ps
ref_p           = 1.0           ; reference pressure, in bar
compressibility = 4.5e-5        ; isothermal compressibility of water, bar^-1
; Add also refcoord-scaling to work with position restraints and pressure coupling
refcoord-scaling = all
; Dispersion correction
DispCorr        = EnerPres      ; account for cut-off vdW scheme
"""

temp_mdp = """
; Temperature coupling is on
tcoupl          = V-rescale     ; modified Berendsen thermostat
tc-grps         = Protein Non-Protein   ; two coupling groups - more accurate
tau_t           = 0.1   0.1     ; time constant, in ps
ref_t           = %s  %s   ; reference temperature, one for each group, in K
"""

vel_mdp = """
; Velocity generation
gen_vel          = yes       ; assign velocities from Maxwell distribution
gen_temp         = %s     ; temperature for Maxwell distribution
gen_seed         = -1        ; generate a random seed
"""

def make_mdp(parms):
  """
  Using GROMACS, we will always solvate, use periodic 
  boundaries and Particle Ewald Mesh.
  """
  mdp = ""
  if 'n_step_minimization' in parms:
    cp = min_mdp
    mdp = cp % parms
  else:
    cp = dynamics_mdp
    mdp += cp % parms
    if 'temp_thermometer' in parms:
      t = parms['temp_thermometer']
      mdp += temp_mdp % (t, t)
    else:
      mdp += "; Temperature coupling is off\n"
      mdp += "tcoupl          = no\n"
    if 'temp_initial' in parms and parms['temp_initial'] > 0.0:
      t = parms['temp_initial']
      mdp += vel_mdp % t
    else:        
      mdp += "; Velocity generation\n"
      mdp += "gen_vel         = no            ; Velocity generation is off \n"
    mdp = "title           = Template for constant temperature/pressure\n" +\
          mdp
  if 'restraint_pdb' in parms and parms['restraint_pdb']:
    mdp += "define            = -DPOSRES  ; position restrain the protein\n"
  return mdp
  
  
def replace_include_file(chain, r1, r2):
  lines = open(chain, 'r').readlines()
  new_lines = []
  for l in lines:
    if l.startswith('#include'):
      new_lines.append(l.replace(r1, r2))
    else:
      new_lines.append(l)
  open(chain, 'w').write(''.join(new_lines))

restraint_header = """
; In this topology include file, you will find position restraint
; entries for all the heavy atoms in your original pdb file.
; This means that all the protons which were added by pdb2gmx are
; not restrained.

[ position_restraints ]
; atom  type      fx      fy      fz
"""

def run(parms):
  name = parms['output_name']
  in_mdp = name + '.in.mdp'
  open(in_mdp, 'w').write(make_mdp(parms))
  mdp = name + '.mdp'

  # Copies across topology files and does
  # some name mangling so that all new topology
  # files point to each other correctly
  top = name + '.top'
  in_top = parms['topology']
  shutil.copy(in_top, top)
  in_name = os.path.basename(in_top).replace('.top', '')
  in_dir = os.path.dirname(in_top)
  file_tag = "%s/%s_*itp" % (in_dir, in_name)
  new_files = [top]
  for f in glob.glob(file_tag):
    new_f = os.path.basename(f)
    new_f = new_f.replace(in_name, name)
    shutil.copy(f, new_f)
    new_files.append(new_f)
  for f in new_files:
    replace_include_file(f, in_name + "_", name + "_")

  restraint_files = [f for f in new_files if 'posre' in f]
  if parms['restraint_pdb']:
    atoms = pdbatoms.Polymer(parms['restraint_pdb']).atoms()
    indices = range(1, len(atoms)+1)
    for fname in restraint_files:
      chain_id = fname.replace('.itp', '').split('_')[-1]
      with open(fname, 'w') as f:
        f.write(restraint_header)
        for i, atom in zip(indices, atoms):
          if atom.chain_id == chain_id and atom.bfactor > 0.0:
            f.write("%6s     1  1000  1000  1000\n" % i)

  in_gro = name + '.in.gro'
  shutil.copy(parms['input_crds'], in_gro)
  tpr = name + '.tpr'
  data.binary(
      'grompp',
      '-f %s -po %s -c %s -p %s -o %s' \
          % (in_mdp, mdp, in_gro, top, tpr),
      tpr)
  util.check_files(tpr)

  data.binary(
      'mdrun',
      '-v -deffnm %s' % (name),
      name + '.mdrun')

  top, crds, vels = get_restart_files(name)
  util.check_files(top, crds)

  delete_backup_files(name)


def merge_simulations(name, pulses):
  save_dir = os.getcwd()
  trr_fname = name + '.trr'
  trrs = [os.path.join(p, trr_fname) for p in pulses]
  trrs = [t for t in trrs if os.path.isfile(t)]
  f = open(trr_fname, 'w')
  for trr in trrs:
    trr = TrrReader(trr)
    for i_frame in range(trr.n_frame-1):
      trr.file.seek(trr.size_frame*i_frame)
      txt = trr.file.read(trr.size_frame)
      f.write(txt)
  f.close()
  for pulse in reversed(pulses):
    gro = "%s/%s.gro" % (pulse, name)
    if os.path.isfile(gro):
      for ext in ['.top', '.itp', '.tpr', '.mdp', '.in.mdp', '.gro']:
        os.system('cp %s/%s*%s .' % (pulse, name, ext))
      os.chdir(save_dir)
      return
  os.chdir(save_dir)
    
    
def preequilibrate(in_name, out_name, temperature):
  top, gro, dummy = get_restart_files(in_name)
  restraint_name = 'restraint'
  parms = langevin_thermometer_parms.copy()
  parms['restraint'] = True
  parms['topology'] = top
  parms['input_crds'] = gro
  parms['output_name'] = restraint_name
  parms['temp_thermometer'] = temperature
  parms['temp_initial'] = temperature
  parms['n_step_per_snapshot'] = 50
  parms['n_step_dynamics'] = 5000
  run(parms)  

  top, gro, dummy = get_restart_files(restraint_name)
  parms = langevin_thermometer_parms.copy()
  parms['topology'] = top
  parms['input_crds'] = gro
  parms['output_name'] = out_name
  parms['temp_thermometer'] = temperature
  parms['n_step_per_snapshot'] = 50
  parms['n_step_dynamics'] = 5000
  run(parms)  


# Trajectory read 


n_dim = 3

class TrrReader:
  def __init__(self, fname):
    self.fname = fname
    self.file = open(self.fname, 'r')

    self.u = xdrlib.Unpacker(self.file.read(200))
    self.magic = self.u.unpack_int()
    self.version = self.u.unpack_string()
    self.size_ir = self.u.unpack_int()
    self.size_e = self.u.unpack_int()
    self.size_box = self.u.unpack_int()
    self.size_vir = self.u.unpack_int()
    self.size_pres = self.u.unpack_int()
    self.size_top = self.u.unpack_int()
    self.size_sym = self.u.unpack_int()
    self.size_x = self.u.unpack_int()
    self.size_v = self.u.unpack_int()
    self.size_f = self.u.unpack_int()
    self.n_atom = self.u.unpack_int()
    self.step = self.u.unpack_int()
    self.nre = self.u.unpack_int()
    self.t = self.u.unpack_float()
    self.lamb = self.u.unpack_float()

    self.calc_precision()
    self.offset = self.u.get_position()
    self.calc_size_frame()

    self.file.seek(0, 2)
    size_file = self.file.tell()
    self.n_frame = size_file / self.size_frame

  def calc_size_frame(self):
    n_vec = 0
    if self.size_box: n_vec += n_dim
    if self.size_vir: n_vec += n_dim
    if self.size_pres: n_vec += n_dim
    if self.size_x: n_vec += self.n_atom
    if self.size_v: n_vec += self.n_atom
    if self.size_f: n_vec += self.n_atom
    self.size_frame = (n_vec*n_dim*self.precision + self.offset)

  def calc_precision(self):
    "Returns 4 for single precision, and 8 for double precision"
    if self.size_box:
      self.precision = self.size_box/n_dim/n_dim
    elif self.size_x:
      self.precision = self.size_x/n_dim
    elif self.size_v:
      self.precision = self.size_v/n_dim
    elif self.size_f:
      self.precision = self.size_f/n_dim
    else:
      raise ValueError, "Cannot determine precision"
    if self.precision not in [4, 8]:
      raise ValueError, "Precision not single or double!"
    
  def next_3_reals(self):
    if self.precision == 4:
      return [self.u.unpack_float() for i in range(3)]
    if self.precision == 8:
      return [self.u.unpack_double() for i in range(3)]

  def __len__(self):
    return self.n_frame
    
  def __repr__(self):
    return "< Gromacs TRR Coord file %s with %d frames of %d atoms >" % \
             (self.trj, self.n_frame, self.n_atom)

  def __getitem__(self, i_frame):
    if i_frame < - 1*self.n_frame or i_frame >= self.n_frame :
      raise IndexError
    elif i_frame < 0 :
      return self.__getitem__(self.n_frame + i_frame)
    self.i_frame = i_frame
    self.file.seek(self.offset + i_frame*self.size_frame)
    self.u = xdrlib.Unpacker(self.file.read(self.size_frame))
    box, positions, velocities, forces = None, None, None, None
    if self.size_box:
      box = [self.next_3_reals() for i in range(n_dim)]
    if self.size_vir:
      dummy = [self.next_3_reals() for i in range(n_dim)]
    if self.size_pres:
      dummy = [self.next_3_reals() for i in range(n_dim)]
    if self.size_x:
      positions = [self.next_3_reals() for i in range(self.n_atom)]
    if self.size_v:
      velocities = [self.next_3_reals() for i in range(self.n_atom)]
    if self.size_f:
      forces = [self.next_3_reals() for i in range(self.n_atom)]
    return box, positions, velocities, forces


class Trajectory(object):
  def __init__(self, name):
    self.name = name
    self.top = name + '.top'
    self.gro = name + '.gro'
    self.trr = name + '.trr'
    self.trr_reader = TrrReader(self.trr)
    self.n_frame = self.trr_reader.n_frame
    self.soup = soup_from_top_gro(self.top, self.gro)
    self.atoms = self.soup.atoms()
    self.load_frame(0)

  def load_frame(self, i_frame):
    if i_frame < -self.n_frame or i_frame >= self.n_frame :
      raise IndexError
    if i_frame < 0 :
      return self.load_frame(self.n_frame + i_frame)
    box, positions, velocities, forces = self.trr_reader[i_frame]
    self.i_frame = i_frame
    n = len(self.atoms)
    for i in range(n):
      v3.set_vector(
          self.atoms[i].pos,
          positions[i][0]*10,
          positions[i][1]*10,
          positions[i][2]*10)
      v3.set_vector(
          self.atoms[i].vel,
          velocities[i][0]*10,
          velocities[i][1]*10,
          velocities[i][2]*10)


