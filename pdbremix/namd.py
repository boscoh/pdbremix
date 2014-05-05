import os
import sys
import copy
import struct
import shutil

import util
import data
import pdbatoms
import v3
import pdbtext

"""
Angstroms/ps in velocity files
"""

# Routines to read and write restart files

def expand_restart_files(name):
  psf = os.path.abspath(name + '.psf')
  coor = os.path.abspath(name + '.coor')
  vel = os.path.abspath(name + '.vel')
  return psf, coor, vel


def get_restart_files(name):
  psf, coor, vel = expand_restart_files(name)
  util.check_files(psf, coor)
  if not os.path.isfile(vel):
    vel = ''
  return psf, coor, vel


def convert_to_namd_atom_names(soup):
  for res in soup.residues():
    if res.type == "ILE" and res.has_atom('CD1'):
      res.change_atom_type('CD1', 'CD')
    if res.has_atom('OXT'):
      res.change_atom_type('OXT', 'OT2')
      if res.has_atom('O'):
        res.change_atom_type('O', 'OT1')
    for atom in res.atoms():
      if atom.type[0].isdigit() and atom.type[1] == "H":
        new_atom_type = atom.type[1:] + atom.type[0]
        res.change_atom_type(atom.type, new_atom_type)
    if res.type == "HOH":
      res.set_type("TIP3")
      res.set_chain_id("W")
    if res.type == "HIS":
      res.set_type("HSE")


def convert_to_pdb_atom_names(soup):
  for res in soup.residues():
    if res.type == "ILE":
      if res.has_atom('CD'):
        res.change_atom_type('CD', 'CD1')
    if res.has_atom('OT2'):
      res.change_atom_type('OT2', 'OXT')
      if res.has_atom('OT1'):
        res.change_atom_type('OT1', 'O')
    if res.type == "TIP3":
      res.set_type("HOH")
      res.set_chain_id(" ")
    if res.type == "HSE":
      res.set_type("HIS")
    for atom in res.atoms():
      if atom.type[-1].isdigit() and atom.type[0] == "H":
        new_atom_type = atom.type[-1] + atom.type[:-1]
        res.change_atom_type(atom.type, new_atom_type)


def read_psf(psf):
  result = []
  reading = False
  for line in open(psf):
    if not reading:
      if "NATOM" in line:
        n_atom = int(line.split()[0])
        reading = True
      continue
    else: # reading
      words = line.split()
      chain_id = words[1]
      q = words[6]
      mass = words[7]
      result.append((chain_id, q, mass))
      if len(result) == n_atom:
        break
  return result


def soup_from_restart_files(psf, in_coor, in_vel):
  soup = pdbatoms.Polymer(in_coor)
  convert_to_pdb_atom_names(soup)
  for atom, (chain_id, q, mass) in zip(soup.atoms(), read_psf(psf)):
    atom.mass = float(mass)
    if chain_id.startswith('WT') or chain_id.startswith('ION'):
      atom.is_hetatm = True
      atom.chain_id = " "
    else:
      atom.chain_id = chain_id[0]
  if in_vel:
    for atom, vel_atom in zip(soup.atoms(), pdbatoms.Polymer(in_vel).atoms()):
      v = vel_atom.pos
      v3.set_vector(atom.vel, v[0], v[1], v[2])
  return soup


def write_soup_to_crds_and_vels(in_soup, name):
  soup = in_soup.copy()
  convert_to_namd_atom_names(soup)
  coor = name + '.coor'
  soup.write_pdb(coor)
  for atom in soup.atoms():
    v3.set_vector(atom.pos, atom.vel[0], atom.vel[1], atom.vel[2])
  vel = name + '.vel'
  soup.write_pdb(vel)
  return coor, vel


def convert_restart_to_pdb(md_name, pdb):
  psf, coor, vel = get_restart_files(md_name)
  soup = soup_from_restart_files(psf, coor, vel)
  soup.write_pdb(pdb)
    


# Routines to generate psf and coor files
  
solvate_script = """
# Solvate
# Set minimum padding
set pad %(solvent_buffer)s

# Run solvate with automatic padding option
package require solvate
resetpsf
solvate %(in_psf)s %(in_pdb)s -o %(name)s.vmd -rotate -rotinc 5 -t $pad

# Find the periodic box size.
mol load psf %(name)s.vmd.psf
set mol1 [molinfo top get id]
animate read pdb %(name)s.vmd.pdb $mol1

set sel1 [atomselect $mol1 all]
set cent1 [measure center $sel1]
set size1 [measure minmax $sel1]
mol delete $mol1

set sizemin [lindex $size1 0]
set sizemax [lindex $size1 1]

set xsize [expr [lindex $sizemax 0]-[lindex $sizemin 0]]
set ysize [expr [lindex $sizemax 1]-[lindex $sizemin 1]]
set zsize [expr [lindex $sizemax 2]-[lindex $sizemin 2]]

# Find the padding in each direction to make the box cubic
set lmax $xsize
if {$ysize > $lmax} {
	set lmax $ysize
}
if {$zsize > $lmax} {
	set lmax $zsize
}

# I like my boxe size to be a nice even number
set maxsize [expr int(ceil($lmax))]
if {$maxsize%%2 != 0} { set maxsize [expr $maxsize +1] }

# Calculate additional padding
set xpad [expr {$pad+0.5*($maxsize-$xsize)}]
set ypad [expr {$pad+0.5*($maxsize-$ysize)}]
set zpad [expr {$pad+0.5*($maxsize-$zsize)}]

puts ":Padding: $xpad $ypad $zpad"
puts ":Box size: $lmax"

# Adjust padding for nonzero center of mass. These are used to manually set the padding in each direction (both + and -)
set xplus [expr $xpad - [lindex $cent1 0]]
set xmin [expr $xpad + [lindex $cent1 0]]
set yplus [expr $ypad - [lindex $cent1 1]]
set ymin [expr $ypad + [lindex $cent1 1]]
set zplus [expr $zpad - [lindex $cent1 2]]
set zmin [expr $zpad + [lindex $cent1 2]]

# Rerun solvate on the original structure using calculated padding to make the box cubic
resetpsf
solvate %(name)s.psf %(name)s.pdb -o %(name)s.vmd -rotate -rotinc 5 -x $xmin +x $xplus -y $ymin +y $yplus -z $zmin +z $zplus

# Check that it worked
mol delete all
mol load psf %(name)s.vmd.psf
set mol1 [molinfo top get id]
animate read pdb %(name)s.vmd.pdb $mol1

set sel1 [atomselect $mol1 all]
set cent1 [measure center $sel1]
set size1 [measure minmax $sel1]

set sizemin [lindex $size1 0]
set sizemax [lindex $size1 1]

set xsize [expr [lindex $sizemax 0]-[lindex $sizemin 0]]
set ysize [expr [lindex $sizemax 1]-[lindex $sizemin 1]]
set zsize [expr [lindex $sizemax 2]-[lindex $sizemin 2]]

puts ":Final size: $xsize, $ysize, $zsize"
puts ":Length: $sizemax"
puts ":Center: $cent1"
puts ":Size: $size1"

# Ionize
package require autoionize
# Find original charge
set all [atomselect $mol1 all]
set q [measure sumweights $all weight charge]

# Determine the number and type of ions to use
set natom [expr round($q)]
if {$natom < 0} {
	set natom [expr abs($natom)]
	autoionize -psf %(name)s.vmd.psf -pdb %(name)s.vmd.pdb -nna $natom -ncl 0 -o %(name)s
} elseif {$natom > 0} {
	autoionize -psf %(name)s.vmd.psf -pdb %(name)s.vmd.pdb -nna 0 -ncl $natom -o %(name)s
} elseif {$natom == 0} {
	exec cp %(name)s.vmd.psf %(name)s.psf
	exec cp %(name)s.vmd.pdb %(name)s.pdb
}

# Check that it worked
mol delete all
mol load psf %(name)s.psf
set mol1 [molinfo top get id]
animate read pdb %(name)s.pdb $mol1
set all [atomselect $mol1 all]
puts ":Old charge: $q, New charge: [measure sumweights $all weight charge]"
mol delete all
"""

def solvate_psf(in_psf, in_pdb, out_name, solvent_buffer=10.0):
  parms = {
    'in_psf': in_psf,
    'in_pdb': in_pdb,
    'name': out_name,
    'solvent_buffer': solvent_buffer,
  }
  tcl = out_name + '.vmd.tcl'
  open(tcl, 'w').write(solvate_script % parms)
  data.binary('vmd', '-dispdev text -eofexit', out_name+'.vmd', tcl)
  util.check_output(out_name + '.vmd.pdb')
  util.check_output(out_name + '.pdb')
 

module_load_script = """
package require psfgen 
topology %(topology)s
alias residue HIS HSE 
alias atom ILE CD1 CD 
"""

fixed_waters_script = """
pdbalias residue HOH TIP3
pdbalias residue WAT TIP3
segment wat { 
 auto none
 pdb %(water_pdb)s
} 
pdbalias atom HOH O OH2 
pdbalias atom WAT O OH2 
coordpdb %(water_pdb)s wat
"""

protein_chain_script = """
segment %(chain_id)s { 
 pdb %(chain_pdb)s 
} 
coordpdb %(chain_pdb)s %(chain_id)s
"""

write_script = """
guesscoord 
writepdb %(out_pdb)s
writepsf %(out_psf)s
"""


def make_disulfide_script(soup):
  script = ""
  n = len(soup.residues())
  for i in range(n):
    for j in range(i+1, n):
      if soup.residue(i).type in 'CYS' and soup.residue(j).type in 'CYS':
        sg1 = soup.residue(i).atom('SG')
        sg2 = soup.residue(j).atom('SG')
        if v3.distance(sg1.pos, sg2.pos) < 3.0:
          chain_id1 = sg1.chain_id
          if chain_id1 == ' ':
            chain_id1 = "A"
          chain_id2 = sg2.chain_id
          if chain_id2 == ' ':
            chain_id2 = "A"
          script += "patch DISU %s:%s %s:%s\n" % (
            chain_id1, i+1, chain_id2, j+1)
  if script:
     script = "# disulfide bonds\n" + script + "\n"
  return script


def pdb_to_top_and_crds(
    force_field, pdb, name,
    solvent_buffer=10.0): 
  """
  Creates charmm .pdb and .psf file. Can only make NAMD
  generate CHARMM topology files, not OPLS ones - that
  requires XPLOR but I don't care. Still, if OPLS
  topology files are provided - can still run.
  """
  solv_dir = name + '.solvate'
  save_dir = os.getcwd()
  util.goto_dir(solv_dir)

  psfgen_psf = name + '.psfgen.psf'
  psfgen_pdb = name + '.psfgen.pdb'
  topology = os.path.join(data.data_dir, 'charmm22.topology')
  stripped_pdb = name + '.clean.pdb' 
  pdbtext.clean_pdb(pdb, stripped_pdb)

  script = module_load_script % { 'topology':topology }
  soup = pdbatoms.Polymer(stripped_pdb)
  for chain_id in soup.chain_ids():
    if chain_id == ' ':
      chain_id = 'A'
    chain_pdb = '%s.chain.%s.pdb' % (name, chain_id) 
    chain = soup.extract_chain(chain_id).write_pdb(chain_pdb)
    script += protein_chain_script % { 
      'chain_id': chain_id, 'chain_pdb': chain_pdb }
  script += make_disulfide_script(soup)
  script += write_script % { 
    'out_pdb': psfgen_pdb, 'out_psf': psfgen_psf }

  psfgen_in = name + ".psfgen.in"
  open(psfgen_in, "w").write(script)
  data.binary('psfgen', psfgen_in, name+'.psfgen')
  util.check_output(psfgen_psf)
  util.check_output(psfgen_pdb)

  solvate_psf(psfgen_psf, psfgen_pdb, name, solvent_buffer)

  psf = name+'.psf'
  coor = name+'.coor'
  pdb = name+'.pdb'
  os.rename(pdb, coor)
  convert_restart_to_pdb(name, pdb)
  shutil.copy(psf, save_dir)
  shutil.copy(coor, save_dir)
  shutil.copy(pdb, save_dir)
  
  os.chdir(save_dir)

  return psf, coor


# Routines to prepare namd minimization/molecular-dynamics runs

minimization_parms = { 
  'topology' : 'in.psf', 
  'input_crds' : 'in.pdb', 
  'output_name' : 'min', 
  'force_field': 'NAMD', 
  'restraint_pdb': '',
  'n_step_minimization' : 100, 
} 

constant_energy_parms = { 
  'topology' : 'in.top', 
  'input_crds': 'in.coor',
  'input_vels': 'in.vel',
  'output_name' : 'md', 
  'random_seed' : 2342, 
  'force_field': 'NAMD', 
  'n_step_per_snapshot' : 5, 
  'n_step_dynamics' : 1000, 
  'restraint_pdb': '',
} 

langevin_thermometer_parms = { 
  'topology' : 'in.top', 
  'input_crds': 'in.coor',
  # use either vel or temp_background to set initial velocities
  'input_vels': '',
  'temp_thermometer' : 300.0, 
  'force_field': 'NAMD', 
  'output_name' : 'md', 
  'random_seed' : 2342, 
  'n_step_per_snapshot' : 5, 
  'n_step_dynamics' : 1000, 
  'restraint_pdb': '',
} 


io_script = """
# load topology and initial coordinates
structure %(topology)s
coordinates %(input_crds)s

# set output basenames
outputName %(output_name)s
binaryOutput no
"""

simulation_parameters_script = """
# forcefield parameters
%(psf_type)s
parameters %(parameter)s

# switching parameters
exclude	scaled1-4
1-4scaling 1.0
dielectric 1.0
switching on
switchDist 8.0
cutoff 12.0
pairListDist 15.0
molly off

# Bond constraints using SHAKE
rigidBonds water

# Periodic Boundary Conditions and PME Electrostatics
wrapWater           on
wrapAll             on
wrapNearest         on
PME                 yes
PMEGridSpacing      1.0
"""

restraint_script="""
# restraints, called constraints here
constraints on
consExp 2
consKFile %(restraint_pdb)s
consKCol B
consRef %(restraint_pdb)s
constraintScaling 100
"""

molecular_dynamics_script = """
# DCD output
dcdFile %(output_name)s.dcd
dcdFreq %(n_step_per_snapshot)s
velDcdFile %(output_name)s.vel.dcd
velDcdFreq %(n_step_per_snapshot)s

# Molecular Dynamics Timesteps
numSteps %(n_step_dynamics)s
timeStep 1  # in fs
firstTimeStep 0
stepsPerCycle 20
nonBondedFreq 1
fullElectFrequency 2
"""

extended_periodic_box_script = """    
# Periodic Box from file
extendedSystem %(xsc)s
"""

new_periodic_box_script = """
# Periodic Box Definition
cellBasisVector1    %(len_x)s     0.0           0.0
cellBasisVector2    0.0           %(len_y)s     0.0
cellBasisVector3    0.0           0.0           %(len_z)s
cellOrigin          %(x_origin)s  %(y_origin)s  %(z_origin)s
"""

def calculate_periodic_box_script(parms):
  script = new_periodic_box_script
  p = pdbatoms.Polymer(parms['input_crds'])
  atoms = p.atoms()
  parms = {}
  for i_v, v in enumerate(['x', 'y', 'z']):
    vals = [a.pos[i_v] for a in atoms]
    v_min, v_max = min(vals), max(vals)
    parms["len_"+v] = v_max - v_min + 0.5
    parms[v+"_origin"] = sum(vals)/float(len(vals))
  return script % parms

import_velocities_script = """
# Import velocities
velocities %(input_vels)s
"""

generate_velocities_script = """
# Generate velocities from temperature
temperature %(temp_initial)s
seed %(random_seed)s
"""

temperature_coupling_script = """
# Temperature Coupling with Langevin Thermometer
langevin on
langevinHydrogen off
langevinTemp %(temp_thermometer)s
langevinDamping 5.0
"""

constant_pressure_script = """
# Constant Pressure Control (variable volume)
useGroupPressure      yes   # yes needed for rigidBonds
useFlexibleCell       no
useConstantArea       no

# Nose-Hoover Langevin barometer
langevinPiston        on
langevinPistonTarget  1.01325   #  in bar -> 1 atm
langevinPistonPeriod  200.0
langevinPistonDecay   100.0
langevinPistonTemp    %(temp_thermometer)s
"""
  
minimization_script = """
# Minimization Parameters
minimization on
numsteps %(n_step_minimization)s
"""

def make_namd_input_file(parms):
  name = parms['output_name']
  script = io_script % parms
  script += simulation_parameters_script % parms
  if parms['restraint_pdb']:
    script += restraint_script % parms
  if parms['xsc']:
    script += extended_periodic_box_script % parms
  else:
    script += calculate_periodic_box_script(parms)
  if 'n_step_dynamics' in parms:
    script += molecular_dynamics_script % parms
    if parms['input_vels']:
      script += import_velocities_script % parms
    elif 'temp_initial' in parms and parms['temp_initial'] > 0.0:
      script += generate_velocities_script % parms
    else:
      raise IndexError("No initial velocity information for dynamics run")
    if 'temp_thermometer' in parms:
      script += temperature_coupling_script % parms
      script += constant_pressure_script % parms
  else:
    script += minimization_script % parms
  return script


def run(in_parms):
  """
  Read parms and creates the appropriate NAMD input files
  for simulation
  """
  parms = copy.deepcopy(in_parms)
  name = parms['output_name']
  namd_in = name + ".namd2.in"

  parms['parameter'] = os.path.join(data.data_dir, 'charmm22.parameter')
  parms['psf_type'] =  'paraTypeCharmm on'

  xsc = parms['topology'].replace('.psf', '.xsc')
  if os.path.isfile(xsc):
    shutil.copy(xsc, name + '.in.xsc')
    parms['xsc'] = name + '.in.xsc'
  else:
    parms['xsc'] = ''
  shutil.copy(parms['topology'], name + '.psf')
  parms['topology'] = name + '.psf'

  shutil.copy(parms['input_crds'], name + '.in.coor')
  parms['input_crds'] = name + '.in.coor'

  if 'input_vels' in parms and parms['input_vels']:
    shutil.copy(parms['input_vels'], name + '.in.vel')
    parms['input_vels'] = name + '.in.vel'
  else:
    parms['input_vels'] = ''

  if 'restraint_pdb' in parms and parms['restraint_pdb']:
    shutil.copy(parms['restraint_pdb'], name + '.restraint.coor')
    parms['restraint_pdb'] = name + '.restraint.coor'
  else:
    parms['restraint_pdb'] = ''
    
  open(namd_in, "w").write(make_namd_input_file(parms))
  
  data.binary('namd2', namd_in, name + '.namd2')

  top, crds, vels = get_restart_files(name)
  util.check_output(top)
  util.check_output(crds)


def merge_simulations(name, pulses):
  for ext in ['.psf', '.coor', '.vel']:
    fname = '%s%s' % (name, ext)
    shutil.copy('%s/%s' % (pulses[-1], fname), fname)
  trajs = [os.path.join(pulse, name + '.dcd') for pulse in pulses]
  merge_trajectories(name + '.psf', trajs, name + '.dcd')
  vels = [os.path.join(pulse, name + '.vel.dcd') for pulse in pulses]
  merge_trajectories(name + '.psf', vels, name + '.vel.dcd')
  

def check_dcd_byte_order(dcd):
  if sys.byteorder in 'big':
    option = '-B'
  elif sys.byteorder in 'little':
    option = '-L'
  else:
    raise "Couldn't find o.s. byte order %s" % sys.byteorder
  data.binary('flipdcd', '%s %s' % (option, dcd), dcd+'.flipdcd')


# Routines to handle the CHARMM dcd trajectories

class DcdReader:
  """
  Read frames from a DCD file in terms of a list (3xn) of floats.

  Data: 
    fname, n_fixed_atom, remarks, 
    free_atom_indices, pos_first_frame, 
    size_frame, n_frame, n_atom

  Methods:
    d = DCD(fname) - open fname and read header
    len(d) - return number of frames
    d[10] - return a frame as a list (3xn) of floats.
  """

  def __init__(self, fname):
    self.fname = fname

    check_dcd_byte_order(self.fname)    
    
    self._file = open(fname, 'rb')

    # Read header, should be 84 bytes
    if self._read_fmt_val('i') != 84:
      raise Exception("DCD: 1st integer is not 84, the size of header")

    if self._read_fmt_vals('4c') != ('C', 'O', 'R', 'D') :
     raise Exception("DCD: Missing CORD tag string")

    self.n_frame = self._read_fmt_val('i') 
    self.i_start = self._read_fmt_val('i')
    self.n_step_save = self._read_fmt_val('i')
    # skip some 
    self._read_fmt_val('5i')
    self.n_fixed_atom = self._read_fmt_val('i')
    self.timeStep = self._read_fmt_val('d')

    self._read_fmt_vals('9i')

    if self._read_fmt_val('i') != 84 :
      raise Exeption("DCD: couldn't find ending 84 of the header")

    # read title block
    size = self._read_fmt_val('i')
    if (size - 4) % 80 != 0 :
      raise "DCD format error 3"
    self.remarks = []
    n_line = self._read_fmt_val('i')
    for i in range(0, n_line):
      s = "".join(self._read_fmt_vals('80c'))
      self.remarks.append(s.strip())
    if self._read_fmt_val('i') != size :
      raise Exception("DCD: not matching title size block end")

    # block to record number of atoms
    if self._read_fmt_val('i') != 4 :
      raise "DCD format error 5"
    self.n_atom = self._read_fmt_val('i')
    if self._read_fmt_val('i') != 4 :
      raise "DCD format error 6"

    if self.n_fixed_atom > 0:
      fmt = "%di" % (self.n_atom - self.n_fixed_atom)
      size = struct.calcsize(fmt)
      if self._read_fmt_val('i') != size:
        raise "DCD format error 7"
      self.free_atom_indices = self._read_fmt_val(fmt)
      for i in range(len(self.free_atom_indices)):
        self.free_atom_indices[i] -= 1
      if self._read_fmt_val('i') != size :
        raise "DCD format error 8"
    else:
      self.free_atom_indices = ()

    self.pos_after_header = self._file.tell()
    size_first_frame = struct.calcsize('%df6i' % (3 * self.n_atom))
    self.pos_first_frame = self._file.tell() + size_first_frame

    n_free_atom = self.n_atom - self.n_fixed_atom
    self.size_frame = struct.calcsize('%df6i' % (3 * n_free_atom))

    if self.n_fixed_atom > 0:
      self._first_frame_x_vals = [0.0 for i in range(self.n_atom)]
      self._first_frame_y_vals = [0.0 for i in range(self.n_atom)]
      self._first_frame_z_vals = [0.0 for i in range(self.n_atom)]
      crds_fmt = '%df' % self.n_atom
      size = struct.calcsize(crds_fmt)

      if size != self._read_fmt_val('i'):
        raise "DCD format error 9"
      self._first_frame_x_vals = self._read_fmt_vals(crds_fmt)
      if size != self._read_fmt_val('i'):
        raise "DCD format error 10"

      if size != self._read_fmt_val('i'):
        raise "DCD format error 11"
      self._first_frame_y_vals = self._read_fmt_vals(crds_fmt)
      if size != self._read_fmt_val('i'):
        raise "DCD format error 12"

      if size != self._read_fmt_val('i'):
        raise "DCD format error 13"
      self._first_frame_z_vals = self._read_fmt_vals(crds_fmt)
      if size != self._read_fmt_val('i'):
        raise "DCD format error 14"

    self._file.seek(0, 2)
    last_pos = self._file.tell()
    size_rest_of_file = last_pos - self.pos_after_header
    implied_size_frame = size_rest_of_file / self.n_frame
    self.extra_block_size = implied_size_frame - self.size_frame 
    if self.extra_block_size > 0:
      self.size_frame += self.extra_block_size

  def _read_fmt_vals(self, fmt):
    return struct.unpack(fmt, self._file.read(struct.calcsize(fmt)))

  def _read_fmt_val(self, fmt):
    return self._read_fmt_vals(fmt)[0]
    
  def __getitem__(self, i):
    if i < - 1*self.n_frame or i >= self.n_frame :
      raise IndexError
    if i < 0 :
      return self.__getitem__(self.n_frame + i)

    if self.n_fixed_atom == 0 :
      pos = self.pos_after_header + i*self.size_frame
      pos += self.extra_block_size
      self._file.seek(pos)

      crds_fmt = "%df" % self.n_atom
      size = struct.calcsize(crds_fmt)
      if size != self._read_fmt_val('i'):
        raise Exception("DCD format error 9")
      x_vals = self._read_fmt_vals(crds_fmt)
      if size != self._read_fmt_val('i'):
        raise Exception("DCD format error 10")
      if size != self._read_fmt_val('i'):
        raise Exception("DCD format error 11")
      y_vals = self._read_fmt_vals(crds_fmt)
      if size != self._read_fmt_val('i'):
        raise Exception("DCD format error 12")
      if size != self._read_fmt_val('i'):
        raise Exception("DCD format error 13")
      z_vals = self._read_fmt_vals(crds_fmt)
      if size != self._read_fmt_val('i'):
        raise Exception("DCD format error 14")

      frame = []
      for x, y, z in zip(x_vals, y_vals, z_vals):
        frame.extend((x, y, z))

    else:
      # first frame cached
      if i == 0:
        return copy.copy(self._firstFrame)

      self._file.seek(self.pos_after_header + (i-1)*self.size_frame)

      crds_fmt = '%df' % len(self.free_atom_indices)
      size = struct.calcsize(crds_fmt)
      if size != self._read_fmt_val('i'):
        raise "DCD format error 9"
      free_x_vals = self._read_fmt_vals(crds_fmt)
      if size != self._read_fmt_val('i'):
        raise "DCD format error 10"
      if size != self._read_fmt_val('i'):
        raise "DCD format error 11"
      free_y_vals = self._read_fmt_vals(crds_fmt)
      if size != self._read_fmt_val('i'):
        raise "DCD format error 12"
      if size != self._read_fmt_val('i'):
        raise "DCD format error 13"
      free_z_vals = self._read_fmt_vals(crds_fmt)
      if size != self._read_fmt_val('i'):
        raise "DCD format error 14"

      x_vals = copy.copy(self._first_frame_x_vals)
      y_vals = copy.copy(self._first_frame_y_vals)
      z_vals = copy.copy(self._first_frame_z_vals)
      k = 0
      for index in self.free_atom_indices:
        x_vals[index] = free_x_vals[k]
        y_vals[index] = free_y_vals[k]
        z_vals[index] = free_z_vals[k]
        k += 1

      frame = []
      for x, y, z in zip(x_vals, y_vals, z_vals):
        frame.extend((x, y, z))

    return frame

  def __len__(self):
    return self.n_frame

  def __del__(self):
    self._file.close()

  def __repr__(self):
    return "< DCD %s with %d frames of %d atoms (%d fixed) >" % \
             (self.fname, self.n_frame, self.n_atom, self.n_fixed_atom)


class Trajectory:
  def __init__(self, name):
    self.name = name
    self.topology = name + '.psf'
    coor = name + '.coor'
    vel = name + '.vel'
    dcd = name + '.dcd'    
    self.coor_traj = DcdReader(dcd)
    vel_dcd = name + '.vel.dcd'
    if os.path.isfile(vel_dcd):
      self.vel_traj = DcdReader(vel_dcd)
    else:
      self.vel_traj = None
    self.soup = soup_from_restart_files(self.topology, coor, vel)
    self.n_frame = len(self.coor_traj)
    self.i_frame = 0
    self.load_frame(self.i_frame)
    
  def load_frame(self, i):
    if i < -self.n_frame or i >= self.n_frame :
      raise IndexError
    if i < 0 :
      return self.load_frame(self.n_frame + i)
    crds = self.coor_traj[i]
    for j, a in enumerate(self.soup.atoms()):
      k = 3*j
      v3.set_vector(a.pos, crds[k], crds[k+1], crds[k+2])
    if self.vel_traj is not None:
      vels = self.vel_traj[i]
      for j, a in enumerate(self.soup.atoms()):
        k = 3*j
        v3.set_vector(a.vel, vels[k], vels[k+1], vels[k+2])
    self.i_frame = i


def merge_trajectories(psf, dcds, out_dcd):
  dcd_reader = DcdReader(dcds[0])
  pos_after_header = dcd_reader.pos_after_header  
  size_frame = dcd_reader.size_frame
  del dcd_reader

  shutil.copy(dcds[0], out_dcd)

  merge_dcd_file = open(out_dcd, "ab+")
  for dcd in dcds[1:]:
    dcd_file = open(dcd, "rb")
    dcd_file.seek(-1, 2)
    eof = dcd_file.tell()

    dcd_file.seek(pos_after_header)
    while dcd_file.tell() < eof:
      merge_dcd_file.write(dcd_file.read(size_frame)) 

    dcd_file.close()
  merge_dcd_file.close()


def preequilibrate(in_name, out_name, temperature):
  top, coor, vel = get_restart_files(in_name)
  restraint_name = 'restraint'
  parms = langevin_thermometer_parms.copy()
  parms['restraint'] = True
  parms['topology'] = top
  parms['input_crds'] = coor
  parms['input_vels'] = vel
  parms['output_name'] = restraint_name
  parms['temp_thermometer'] = temperature
  parms['temp_initial'] = temperature
  parms['n_step_per_snapshot'] = 50
  parms['n_step_dynamics'] = 1000
  run(parms)  

  top, coor, vel = get_restart_files(restraint_name)
  parms = langevin_thermometer_parms.copy()
  parms['topology'] = top
  parms['input_crds'] = coor
  parms['input_vels'] = vel
  parms['output_name'] = out_name
  parms['temp_thermometer'] = temperature
  parms['n_step_per_snapshot'] = 50
  parms['n_step_dynamics'] = 1000
  run(parms)  

