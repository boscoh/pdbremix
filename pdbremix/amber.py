"""
This modules interaces with the amber molecular dynamics package,
which should already be installed and available on the path. The
routines are used to access the topology files and get the frames
of a trajectory. The interface is mainly through pdb files and
strings.

valid force-fields are AMBER11 and AMBER8
"""


import os
import copy
import shutil
import string
import re

import util
import v3
import pdbtext
import pdbatoms
import data


# Routines to handle topology coordinate, and restart files

def expand_restart_files(name):
  top = os.path.abspath(name + '.top')
  crds = os.path.abspath(name + '.rst')
  if not os.path.isfile(crds):
    crds = os.path.abspath(name + '.crd')
  vels = ''
  return top, crds, vels


def get_restart_files(name):
  top, crds, vels = expand_restart_files(name)
  util.check_files(top, crds)
  return top, crds, vels


def read_top(top):
  section = None
  len_field = None
  parse = None
  parse_map = { 'a':str, 'I':int, 'E':float }

  topology = {}
  for line in open(top, "rU"):
    line = line[:-1]
    if line.startswith("%"):
      words = line.split()
      key = words[0][1:]
      if key == "FLAG":
        section = words[1]
        topology[section] = []
      elif key.startswith("FORMAT"):
        format_str = key[7:-1]
        len_field = int(re.split(r'\D+', format_str)[1])
        val_type = re.search('(a|I|E)', format_str).group(0)
        parse = parse_map[val_type]
    else:
      indices = range(0, len(line), len_field)
      pieces = [line[i:i+len_field] for i in indices]
      topology[section].extend(map(parse, pieces))

  name_str = """
  NATOM NTYPES NBONH MBONA NTHETH 
  MTHETA NPHIH MPHIA NHPARM NPARM 
  NNB NRES NBONA NTHETA NPHIA 
  NUMBND NUMANG NPTRA NATYP NPHB 
  IFPERT NBPER NGPER NDPER MBPER 
  MGPER MDPER IFBOX NMXRS IFCAP
  NUMEXTRA NCOPY  """
  for name, val in zip(name_str.split(), topology['POINTERS']):
    topology[name] = val

  return topology


def soup_from_topology(topology):
  soup = pdbatoms.Polymer()
  chain_id = ''
  n_res = topology['NRES']
  n_atom = topology['NATOM']
  for i_res in range(n_res):
    res_type = topology['RESIDUE_LABEL'][i_res].strip()
    if res_type == "WAT":
      res_type = "HOH"
    res = pdbatoms.Residue(res_type, chain_id, i_res+1)
    soup.append_residue(res)
    res = soup.residue(i_res)
    i_atom_start = topology['RESIDUE_POINTER'][i_res] - 1
    if i_res == n_res-1:
      i_atom_end = n_atom
    else:
      i_atom_end = topology['RESIDUE_POINTER'][i_res+1] - 1
    for i_atom in range(i_atom_start, i_atom_end):
      atom = pdbatoms.Atom()
      atom.vel = v3.vector()
      atom.num = i_atom+1
      atom.res_num = i_res+1
      atom.res_type = res_type
      if atom.res_type == data.solvent_res_types:
        print 'hetatm'
        atom.is_hetatm = True
      atom.type = topology['ATOM_NAME'][i_atom].strip()
      atom.mass = topology['MASS'][i_atom]
      atom.charge = topology['CHARGE'][i_atom]/18.2223
      atom.element = pdbatoms.guess_element(
          atom.res_type, atom.type)
      soup.insert_atom(-1, atom)
  # TODO detect chains
  return soup


def convert_to_pdb_atom_names(soup):
  for res in soup.residues():
    if res.type in data.solvent_res_types:
      for a in res.atoms():
        a.is_hetatm = True
    if res.type == "HSE":
      res.set_type("HIS")
    for atom in res.atoms():
      if atom.type[-1].isdigit() and atom.type[0] == "H":
        new_atom_type = atom.type[-1] + atom.type[:-1]
        res.change_atom_type(atom.type, new_atom_type)


def load_crd_or_rst_into_soup(soup, crd_or_rst):
  f = open(crd_or_rst, "r")
  f.readline()

  n_atom = int(f.readline().split()[0])
  n_crd = n_atom * 3
  n_line = n_crd / 6
  if n_crd % 6 > 0:
    n_line += 1

  line_list = [f.readline()[:-1] for i in range(0, n_line)]
  s = "".join(line_list)
  vals = [float(s[i:i+12]) for i in xrange(0, len(s), 12)]
  if len(vals) != n_crd:
    raise ValueError, "Improper number of coordinates in rst file."
    
  for i, atom in enumerate(sorted(soup.atoms(), pdbatoms.cmp_atom)):
    v3.set_vector(atom.pos, vals[i*3], vals[i*3+1], vals[i*3+2])

  if crd_or_rst.endswith('.rst'):
    line_list = [f.readline()[:-1] for i in range(0, n_line)]
    s = "".join(line_list)
    vals = [float(s[i:i+12]) for i in xrange(0, len(s), 12)]
    if len(vals) != n_crd:
      raise ValueError, "Improper number of coordinates in rst file."

    convert_amber_vel_to_angstroms_per_ps = 20.455
    for i, atom in enumerate(sorted(soup.atoms(), pdbatoms.cmp_atom)):
      v3.set_vector(atom.vel, vals[i*3], vals[i*3+1], vals[i*3+2])
      atom.vel = v3.scale(atom.vel, convert_amber_vel_to_angstroms_per_ps)

  f.close()


def strip_cr(line):
  if line.endswith('\n'):
    return line[:-1]
  if line.endswith('\r\n'):
    return line[:-1]
  return line


def soup_from_top_and_crd_or_rst(top, crd_or_rst):
  topology = read_top(top)
  soup = soup_from_topology(topology)
  convert_to_pdb_atom_names(soup)
  load_crd_or_rst_into_soup(soup, crd_or_rst)
  if topology['IFBOX'] > 0:
    lines = open(crd_or_rst, "r").readlines()
    lines = [l for l in reversed(lines) if l.strip()]
    soup.box_dimension_str = strip_cr(lines[0])
  return soup


def soup_from_restart_files(top, crds, vels):
  return soup_from_top_and_crd_or_rst(top, crds)


def convert_restart_to_pdb(md_name, pdb):
  top, crds, vels = get_restart_files(md_name)
  soup = soup_from_restart_files(top, crds, vels)
  soup.write_pdb(pdb)


def write_soup_to_rst(soup, rst):
  f = open(rst, "w")

  f.write(" ".ljust(80) + "\n")
  f.write("%5d  0.0000000E+00\n" % len(soup.atoms()))

  i = 0
  for atom in sorted(soup.atoms(), pdbatoms.cmp_atom):
    x, y, z = atom.pos
    f.write("%12.7f%12.7f%12.7f" % (x, y, z))
    i += 1
    if i % 2 == 0:
      f.write("\n")
      i = 0
  if len(soup.atoms()) % 2 != 0:
    f.write("\n")

  i = 0
  convert_to_amber_vel = 1.0 / 20.455
  for atom in sorted(soup.atoms(), pdbatoms.cmp_atom):
    x, y, z = atom.vel
    vx = x * convert_to_amber_vel
    vy = y * convert_to_amber_vel
    vz = z * convert_to_amber_vel
    f.write("%12.7f%12.7f%12.7f" % (vx, vy, vz))
    i += 1
    if i % 2 == 0:
      f.write("\n")
  if len(soup.atoms()) % 2 != 0:
    f.write("\n")

  if hasattr(soup, 'box_dimension_str') and soup.box_dimension_str is not None:
    f.write(soup.box_dimension_str + "\n")
    
  f.close()
  
  
def write_soup_to_crds_and_vels(soup, name):
  write_soup_to_rst(soup, name + '.rst')
  return name + '.rst', ''
  

# Routines to genearte topology and coordinate files from PDB

force_field_script = """
# leaprc to generate AMBER topology and coordinate

# load in amber force field
source %(amber_ff)s

# use AMBER6 PB radii as we will use igb=1, gbparm=2
set default PBradii mbondi
"""

explicit_water_box_script = """
# add explicit waters
solvateBox pdb TIP3PBOX %(solvent_buffer)f iso
"""

save_and_quit_script = """
# save topology and coordinates
saveAmberParm pdb %(top)s %(crd)s
quit
"""

def disulfide_script_and_rename_cysteines(in_pdb, out_pdb):
  soup = pdbatoms.Polymer(in_pdb)
  script = " # disulfide bonds\n"
  n = len(soup.residues())
  for i in range(n):
    for j in range(i+1, n):
      if soup.residue(i).type in 'CYS' and soup.residue(j).type in 'CYS':
        p1 = soup.residue(i).atom('SG').pos
        p2 = soup.residue(j).atom('SG').pos
        if v3.distance(p1, p2) < 3.0:
          soup.residue(i).type = 'CYX'
          for atom in soup.residue(i).atoms():
            atom.res_type = 'CYX'
          soup.residue(j).type = 'CYX'
          for atom in soup.residue(j).atoms():
            atom.res_type = 'CYX'
          script += "bond pdb.%d.SG pdb.%d.SG\n" % (i+1, j+1)
  soup.write_pdb(out_pdb)
  return script


def run_tleap(
    force_field, pdb, name, solvent_buffer=0.0, excess_charge=0): 
  "Convert a .pdb file into amber .top and .crd file."
  util.check_output(pdb)
  tleap_pdb = name + '.clean.pdb'
  pdbtext.clean_pdb(pdb, tleap_pdb)
  util.check_output(tleap_pdb)

  top = name + '.top'
  crd = name + '.crd'

  params = { 
    'top': top, 
    'crd': crd, 
    'pdb': tleap_pdb,
    'data_dir':data.data_dir,
    'solvent_buffer': solvent_buffer,
  }

  if 'AMBER11' in force_field:
    params['amber_ff'] = "leaprc.ff99SB"
  elif 'AMBER8' in force_field:
    params['amber_ff'] = "leaprc.ff96"
  else:
    raise Exception("Don't know which version of AMBER(8|11) to use.")

  script = force_field_script
  residues = [r.type for r in pdbatoms.Polymer(tleap_pdb).residues()]
  if 'PHD' in residues:
    leaprc = open("%s/phd.leaprc" % data.data_dir).read()
    script += leaprc
  if 'ZNB' in residues:
    leaprc = open("%s/znb.leaprc" % data.data_dir).read()
    script += leaprc
  script += "pdb = loadpdb %(pdb)s\n"
  script += disulfide_script_and_rename_cysteines(tleap_pdb, tleap_pdb)
  util.check_output(tleap_pdb)

  if 'GBSA' not in force_field:
    if excess_charge != 0:
      if excess_charge > 0:
        script += "addions pdb Cl- 0\n"
      else:
        script += "addions pdb Na+ 0\n"
    solvent_buffer = 10
    params['solvent_buffer'] = solvent_buffer
    script += explicit_water_box_script

  script += save_and_quit_script

  script = script % params

  tleap_in = name + ".tleap.in"
  tleap = data.binary("tleap")
  open(tleap_in, "w").write(script)
  data.binary('tleap', "-f "+tleap_in, name+'.tleap')
  if os.path.isfile('leap.log'):
    os.rename('leap.log', name + '.tleap.log')
  util.check_output(name+'.tleap.log', ['FATAL'])
  util.check_output(top)
  util.check_output(crd)

  return top, crd
  

def pdb_to_top_and_crds(
    force_field, pdb, name, solvent_buffer=0.0): 
  top, crd = run_tleap(
      force_field, pdb, name, solvent_buffer, 0)
  if 'GBSA' not in force_field:
    # redo for counterions
    charges = read_top(name+'.top')['CHARGE']
    charge = int(round(sum(charges)/18.2223))
    if charge != 0:
      top, crd = run_tleap(
          force_field, pdb, name, solvent_buffer, charge)

  convert_restart_to_pdb(name, name+'.pdb')

  return top, crd


# Simulation routines

minimization_parms = { 
  'topology' : 'in.top', 
  'input_crds' : 'in.crd', 
  'output_name' : 'min', 
  'force_field': 'GBSA',
  'restraint_pdb': '',
  'n_step_minimization' : 100, 
} 

constant_energy_parms = { 
  'topology' : 'in.top', 
  'input_crds' : 'in.crd', 
  'output_name' : 'md', 
  'force_field': 'GBSA',
  'restraint_pdb': '',
  'n_step_per_snapshot' : 5, 
  'n_step_dynamics' : 1000, 
} 

langevin_thermometer_parms = { 
  'topology' : 'in.top', 
  'input_crds' : 'in.crd', 
  'output_name' : 'md', 
  'force_field': 'GBSA',
  'restraint_pdb': '',
  'random_seed' : 2342, 
  'temp_thermometer' : 300.0, 
  'temp_initial': 0.0, # ignored if it is 0.0
  'n_step_per_snapshot' : 5, 
  'n_step_dynamics' : 1000, 
  'n_step_per_thermostat' : 100, 
} 


# frequent low energy calculation
sander_script = """
generated by amber.py
&cntrl
"""


# no periodicity, generatlized born, and surface area terms
gbsa_script = "  ntb = 0, igb = 2, gbsa = 1, cut = 12.0,"


# peridoicity/constant pressure, isotropic position scaling, no gb/sa
explicit_water_script = "  ntb = 2, ntp = 1, igb = 0, gbsa = 0, cut = 8.0,"


# 10 steps of steepest descent then conjugate gradient for rest of steps
minimization_script = """
  imin = 1, ntmin = 1, maxcyc = %(n_step_minimization)s, ncyc = 10,
"""


molecular_dynamics_script = """
  ntpr = %(n_step_per_snapshot)s, ntave = 0, ntwr = 500, iwrap = 0, ioutfm = 0,
  ntwx = %(n_step_per_snapshot)s, ntwv = %(n_step_per_snapshot)s, 
  ntwe = %(n_step_per_snapshot)s, 
  nstlim = %(n_step_dynamics)s, nscm = 50, dt = 0.001,
  nrespa = 1,
"""

# langevin thermometer
thermostat_script = """
  ntt = 3, gamma_ln = 5, temp0 = %(temp_thermometer)s, vlimit = 0.0,
  ig = %(random_seed)s, 
  tempi = %(temp_initial)s, 
"""


def make_sander_input_file(parms):
  soup = soup_from_top_and_crd_or_rst(
      parms['topology'], parms['input_crds'])

  has_periodic_box = \
      hasattr(soup, 'box_dimension_str') and \
      soup.box_dimension_str

  script = sander_script % parms

  # all bonds calculated
  script += "  ntf = 1, ntc = 1,\n"

  if parms['restraint_pdb']:
    script += "  ntr = 1,\n"

  if has_periodic_box:
    script += explicit_water_script
  else:
    script += gbsa_script

  if 'n_step_minimization' in parms:
    script += minimization_script 
  elif 'n_step_dynamics' in parms:
    if parms['input_crds'].endswith('.rst'):
      script += "  ntx = 5, irest = 1,\n"
    else:
      script += "  ntx = 1,\n"
    script += molecular_dynamics_script
    if 'temp_thermometer' in parms:
      script += thermostat_script
  else:
    raise Exception("Can't parse parameters to run")

  script += "&end\n"

  return script % parms


restraint_script = """FIND
* * S *
* * B *
* * 3 *
* * E *
SEARCH
"""


def make_restraint_script(pdb):
  util.check_output(pdb)
  script = "Restrained atoms from %s\n" % pdb
  restraint_weight = 100.0
  script += "%s\n" % restraint_weight
  script += restraint_script
  for i, atom in enumerate(pdbatoms.AtomList(pdb).atoms()):
    if atom.bfactor > 0.0:
      script += "ATOM %d %d\n" % (i+1, i+1)
  script += "END\n"
  script += "END\n"
  return script


def run(in_parms):
  "Runs AMBER with in_parms"
  parms = copy.deepcopy(in_parms)
  name = parms['output_name']

  input_top = parms['topology'] 
  util.check_files(input_top)
  new_top = name + '.top'
  shutil.copy(input_top, new_top)

  input_crd = parms['input_crds']
  util.check_files(input_crd)
  if input_crd.endswith('.crd'): 
    new_crd = name + '.in.crd'
  else:
    new_crd = name + '.in.rst'
  shutil.copy(input_crd, new_crd)
  
  if 'n_step_minimization' in parms:
    rst = name + ".crd"
  else:
    rst = name + ".rst"
  trj = name + ".trj"
  vel_trj = name + ".vel.trj"
  ene = name + ".ene"
  inf = name + ".inf"
  sander_out = name + ".sander.out"
  sander_in = name + ".sander.in"

  sander = data.binary("sander")
  script = make_sander_input_file(parms)
  args = "-O -i %s -o %s -p %s -c %s -r %s -x %s -v %s -e %s -inf %s" \
          % (sander_in, sander_out, new_top, new_crd, rst, trj, vel_trj, ene, inf)

  if parms['restraint_pdb']:
    pdb = parms['restraint_pdb']
    script += make_restraint_script(pdb)
    soup = pdbatoms.Polymer(pdb)
    ref_crd = name + '.restraint.crd'
    write_soup_to_rst(soup, ref_crd)
    util.check_output(ref_crd)
    args += " -ref %s" % ref_crd

  open(sander_in, "w").write(script)

  data.binary('sander', args, name)

  util.check_output(sander_out, ['FATAL'])
  top, crds, vels = get_restart_files(name)
  util.check_output(top)
  util.check_output(crds)


def low_temperature_equilibrate(in_name, out_name, temperature):
  """
  The problem with equilibration is that since at low
  temperatures the simulations are run at constant energy, we
  need to prequilibrate at constant energy as well. With
  implicit solvation, this is an unstable process as the system
  quickly rises to 20 K regardless of the starting conformation.
  We need to let the system find a better minimum *at constant
  energy*.
  """

  top, crd, dummy = get_restart_files(in_name)
  md_dirs = ['heat1', 'const2', 'heat3']

  util.goto_dir(md_dirs[0])
  parms = langevin_thermometer_parms.copy()
  parms['topology'] = top
  parms['input_crds'] = crd
  parms['output_name'] = 'md'
  parms['temp_thermometer'] = temperature
  parms['temp_initial'] = temperature
  parms['n_step_per_snapshot'] = 50
  parms['n_step_dynamics'] = 1000
  run(parms)
  in_top, in_crds, in_vels = get_restart_files('md')

  util.goto_dir(os.path.join('..', md_dirs[1]))
  parms = constant_energy_parms.copy()
  parms['topology'] = top
  parms['input_crds'] = crd
  parms['output_name'] = 'md'
  parms['n_step_per_snapshot'] = 50
  parms['n_step_dynamics'] = 10000
  run(parms)
  in_top, in_crds, in_vels = get_restart_files('md')

  util.goto_dir(os.path.join('..', md_dirs[2]))
  parms = langevin_thermometer_parms.copy()
  parms['topology'] = top
  parms['input_crds'] = crd
  parms['output_name'] = 'md'
  parms['temp_thermometer'] = temperature
  parms['temp_initial'] = temperature
  parms['n_step_per_snapshot'] = 50
  parms['n_step_dynamics'] = 1000
  run(parms)
  
  util.goto_dir('..')
  merge_simulations('md', md_dirs)


def calculate_energy(top, crd):
  "Returns potential energy of top and crd by running sander"
  
  top = os.path.abspath(top)
  crd = os.path.abspath(crd)

  util.goto_dir('energy-temp')

  parms = minimization_parms.copy()
  parms['topology'] = top
  parms['input_crds'] = crd
  parms['output_name'] = 'energy'
  parms['n_step_minimization'] = 0
  parms['n_step_steepest_descent'] = 0

  run(parms)

  lines = open('energy.sander.out').readlines()

  energies = {}
  is_results = False
  for i, line in enumerate(lines):
    if not is_results:
      if '4' in line and 'RESULTS' in line:
        is_results = True
    else:
      if 'NSTEP' in line and 'ENERGY' in line:
        util.goto_dir('..')
        util.clean_fname('energy-temp')
        return float(lines[i+1].split()[1])

  raise "Can't find energies"


# Routines to load trajctories

class TrjReader:
  """
  Class that reads successive sets of coordinates from Amber trajectory (or
  velocity) files. Note: for coordinate files, the box dimensions (3 floats)
  are included at the end of each frame, so this is taken care of in the code
  to read the frame, in order to calculate the correct frame size. However,
  the box dimensions are *not* included in velocity coordinate files.
  """
  
  def __init__(self, top, trj):
    """
    trj: name of amber trajectory (or velocity) file containing sets of
         coordinates 
    top: name of amber top file of trajectory
    """
    
    self.trj = trj
    self.top = top

    self.topology = read_top(top)

    if '.trj' in trj:
      self.is_skip_box_dims = self.topology['IFBOX'] > 0
    else:
      self.is_skip_box_dims = False

    if self.trj.split(".")[-1].strip().lower() == "gz":
      self._file_open = gzip.GzipFile
    else:
      self._file_open = open
    
    self.file = self._file_open(self.trj, "r")

    # skip header line
    self.pos_start_frame = len(self.file.readline())

    # calculate the size of each frame
    self.n_atom = self.topology['NATOM']
    if self.n_atom == 0:
      raise "No names found in .top file"

    n_line = (3 * self.n_atom) / 10
    if (3 * self.n_atom) % 10 > 0: 
      n_line += 1

    if self.is_skip_box_dims:
      n_line += 1

    self.size_frame = 0
    for i in range(0, n_line):
      self.size_frame += len(self.file.readline())
      
    # calculate n_frame
    self.file.seek(0, 2)
    pos_eof = self.file.tell()
    self.n_frame = int((pos_eof - self.pos_start_frame) / self.size_frame)

  def __getitem__(self, i):
    "Gets the list of floats for frame i."

    # Find place in file
    if i < - 1*self.n_frame or i >= self.n_frame :
      raise IndexError
    elif i < 0 :
      return self.__getitem__(self.n_frame + i)
    else :
      self.file.seek(self.pos_start_frame + i*(self.size_frame))

    # read frame as a string
    s = self.file.read(self.size_frame).replace("\n", "")

    # convert string to float
    vals = [float(s[i:i+8]) for i in xrange(0, len(s), 8)]
    
    if self.is_skip_box_dims:
      vals = vals[:-3]

    if len(vals) != 3 * self.n_atom:
      raise ValueError, "Improper number of coordinates in frame."
    
    return vals

  def save_to_crd(self, crd, i_frame):
    "Saves coordinates of i_frame to an AMBER crd file."
    coords = self[i_frame]

    f = open(crd, "w")

    f.write("ACE".ljust(80) + "\n")
    f.write("%5d  0.0000000E+00\n" % (len(coords) // 3))

    p = ["%12.7f" % x for x in coords]

    n_val_per_line = 6
    r = len(p) % n_val_per_line
    if r > 0: 
      p.extend([""] * (n_val_per_line - r))

    for i in xrange(0, len(p), n_val_per_line):
      f.write("".join(p[i:i + n_val_per_line]) + "\n")

    f.close()  

  def __len__(self):
    return self.n_frame
    
  def __del__(self):
    self.file.close()

  def __repr__(self):
    return "< Amber Coord file %s with %d frames of %d atoms >" % \
             (self.trj, self.n_frame, self.n_atom)
    

def convert_crd_to_trj_frame(crd):
  vals = [float(word) for word in util.words_in_file(crd)[1:]]
  lines = []
  line = ''
  for i in range(0, len(vals)):
    line += "%8.3f" % vals[i]
    if (i % 10) == 9:
      lines.append(line)
      line = ''
  if line:
    lines.append(line)
  return '\n'.join(lines) + '\n'


class Trajectory:
  """
  Class to hold chains of a pdb that changes according to a dcd
  
  Data:
    soup: collection of chains that get updated
    n_frame: number of frames in trajectory
  
  Methods:
    load_frame: loads the i'th frame from dcdectory
  """
  def __init__(self, name):
    self.name = name
    self.top = name + '.top'

    coor_trj_fname = name + '.trj'
    self.coor_traj = TrjReader(self.top, coor_trj_fname)

    vel_trj_fname = name + '.vel.trj'
    if os.path.isfile(vel_trj_fname):
      self.vel_traj = TrjReader(self.top, vel_trj_fname)
    else:
      self.vel_traj = None

    self.coor_traj.save_to_crd(self.name+'temp.crd', 0)
    self.soup = soup_from_topology(self.coor_traj.topology)
    load_crd_or_rst_into_soup(self.soup, self.name+'temp.crd')
    util.clean_fname(self.name+'temp.crd')

    self.n_frame = len(self.coor_traj)
    self.i_frame = 0
    self.load_frame(self.i_frame)

  def load_frame(self, i):
    if i < -self.n_frame or i >= self.n_frame :
      raise IndexError
    if i < 0 :
      return self.load_frame(self.n_frame + i)
    coords = self.coor_traj[i]
    for j, a in enumerate(self.soup.atoms()):
      k = 3*j
      v3.set_vector(a.pos, coords[k], coords[k+1], coords[k+2])
    if self.vel_traj is not None:
      vels = self.vel_traj[i]
      for j, a in enumerate(self.soup.atoms()):
        k = 3*j
        v3.set_vector(a.vel, vels[k], vels[k+1], vels[k+2])
    self.i_frame = i


def merge_trajectories(top, trajs, out_traj):
  """
  Given a list of traj filenames (trajs), merges them into one complete
  trajectory (out_traj) using top to work out the number of atoms, and
  hence the size of the frame of the trajectory.
  """
  
  trj_reader = TrjReader(top, trajs[0])
  pos_start_frame = trj_reader.pos_start_frame  
  size_frame = trj_reader.size_frame
  del trj_reader

  shutil.copy(trajs[0], out_traj)

  merge_traj_file = open(out_traj, "ab+")
  for traj in trajs[1:]:
    traj_file = open(traj, "rb")
    traj_file.seek(-1, 2)
    eof = traj_file.tell()
    traj_file.seek(pos_start_frame)
    while traj_file.tell() < eof:
      merge_traj_file.write(traj_file.read(size_frame)) 
    traj_file.close()
  merge_traj_file.close()


def merge_simulations(name, pulses):
  """
  Given a list of directories with partial trajectories in each directory
  with the same name for the md, will splice them together into one uber
  simulation.
  """
  
  shutil.copy(os.path.join(pulses[0], name + '.sander.in'), name + '.sander.in')
  shutil.copy(os.path.join(pulses[0], name + '.top'), name + '.top')
  shutil.copy(os.path.join(pulses[-1], name + '.rst'), name + '.rst')

  # merge energies of pulses into one energy file
  f = open(name + '.energy', 'w')
  f.write('[\n')
  n_step = 0
  time = 0.0
  for pulse in pulses:
    energy_fname = os.path.join(pulse, name + '.energy')
    if os.path.isfile(energy_fname):
      blocks = eval(open(energy_fname).read())
    else:
      sander_out = os.path.join(pulse, name + '.sander.out')
      blocks = read_time_blocks(sander_out)
      for block in blocks:
        block_n_step = int(block['NSTEP'])
        block_time = float(block['TIME(PS)'])
        block['NSTEP'] = str(block_n_step + n_step)
        block['TIME(PS)'] = str(block_time + time)
        f.write(str(block) + ',\n')
    n_step = int(blocks[-1]['NSTEP'])
    time = float(blocks[-1]['TIME(PS)'])
  f.write(']\n')
  f.close()
    
  trajs = [os.path.join(pulse, name + '.trj') for pulse in pulses]
  merge_trajectories(name + '.top', trajs, name + '.trj')

  vels = [os.path.join(pulse, name + '.vel.trj') for pulse in pulses]
  merge_trajectories(name + '.top', vels, name + '.vel.trj')
  

def read_time_blocks(fname):
  "Returns a list of dictionaries containing energy values from sander out file"

  allowed_keys = ['EKtot', 'Etot', 'TEMP(K)', 
                  'TIME(PS)', 'EPtot', 'NSTEP']

  def extract_to_dict(line, ext_dict):
    words = line.split()
    for i in range(len(words)):
      if words[i] == "=":
        key, value = words[i-1], words[i+1]
        if key in allowed_keys:
          ext_dict[key] = value

  blocks = []
  block_dict = {}
  is_results = False
  for line in open(fname):
    if is_results:
      if 'A V E R A G E S' in line:
        is_results = False
      elif '----' in line:
        if len(block_dict.keys()) > 0:
          blocks.append(block_dict.copy())
      else:
        extract_to_dict(line, block_dict)    
    else:
      if '4' in line and 'RESULTS' in line:
        is_results = True
        i = -1
  return blocks


