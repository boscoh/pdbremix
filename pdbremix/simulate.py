import os
import shutil
import copy 

import util
import amber
import namd
import gromacs
import data


def get_md_module(force_field):
  if force_field.startswith('GROMACS'):
    return gromacs
  elif force_field.startswith('AMBER'):
    return amber
  elif force_field.startswith('NAMD'):
    return namd
  else:
    raise ValueError, "unrecognized force-field" + force_field


def pdb_to_top_and_crds(
    force_field, raw_pdb, md_name, solvent_buffer=10.0):
  "Returns topology and coordinate fnames"
  md_module = get_md_module(force_field)
  top, crd = md_module.pdb_to_top_and_crds(
      force_field, raw_pdb, md_name, solvent_buffer)
  return os.path.abspath(top), os.path.abspath(crd)
    

def fetch_parms(
    force_field, top, crds, restraint_pdb, 
    simulation_type, md_name):
  """
  Get a dictionary with the necessary fields to run a simulation.
  Parms: simulation_type can be minimization, langevin_thermometer,
         constant_energy
  """
  md_module = get_md_module(force_field)
  parms_dict_name = '%s_parms' % simulation_type
  parms = getattr(md_module, parms_dict_name).copy()
  parms['force_field'] = force_field
  parms['topology'] = top
  parms['input_crds'] = crds
  parms['output_name'] = md_name
  parms['cutoff'] = 12.0
  if restraint_pdb:
    parms['restraint_pdb'] = restraint_pdb
  return parms


def expand_restart_files(force_field, md_name):
  md_module = get_md_module(force_field)
  return md_module.expand_restart_files(md_name)


def get_restart_files(md_name):
  for module in [amber, gromacs, namd]:
    try:
      return module.get_restart_files(md_name)
    except util.FileException:
      pass
  raise Exception("Couldn't find restart files for " + md_name)
  

def convert_restart_to_pdb(md_name, pdb):
  for module in [amber, gromacs, namd]:
    try:
      return module.convert_restart_to_pdb(md_name, pdb)
    except util.FileException:
      pass
  raise Exception("Couldn't find restart files for " + md_name)


def write_soup_to_crds_and_vels(force_field, soup, md_name):
  md_module = get_md_module(force_field)
  return md_module.write_soup_to_crds_and_vels(soup, md_name)


def soup_from_restart_files(top, crds, vels, skip_solvent=True):
  if crds.endswith('.gro'):
    return gromacs.soup_from_restart_files(
        top, crds, vels, skip_solvent)
  elif top.endswith('.top'):
    return amber.soup_from_restart_files(top, crds, vels)
  elif top.endswith('.psf'):
    return namd.soup_from_restart_files(top, crds, vels)


def run_parms(parms):
  md_module = get_md_module(parms['force_field'])
  md_module.run(parms)
  

def minimize(
    force_field, top, crds, md_name, 
    restraint_pdb="", n_step=200):
  parms = fetch_parms(
      force_field, top, crds, restraint_pdb, 
      'minimization', md_name)
  parms['n_step_minimization'] = n_step
  run_parms(parms)


def langevin(
    force_field, top, crds, vels, n_step, temp, md_name, 
    n_step_per_snapshot=50, restraint_pdb=""):
  parms = fetch_parms(
      force_field, top, crds, restraint_pdb, 
      'langevin_thermometer', md_name)
  parms['input_vels'] = vels
  parms['n_step_dynamics'] = n_step
  parms['n_step_per_snapshot'] = n_step_per_snapshot
  parms['temp_thermometer'] = "%.1f" % temp
  parms['temp_initial'] = str(temp)
  parms['gamma_ln'] = 5.0
  run_parms(parms)


def constant(
    force_field, top, crds, vels, n_step, md_name, 
    n_step_per_snapshot=50, restraint_pdb=""):
  parms = fetch_parms(
      force_field, top, crds, restraint_pdb, 'constant_energy', md_name)
  parms['input_vels'] = vels
  parms['temp_thermometer'] = 10.0
  parms['temp_initial'] = 10.0
  parms['n_step_dynamics'] = n_step
  parms['cutoff'] = 12.0
  parms['n_step_per_snapshot'] = n_step_per_snapshot
  run_parms(parms)


def merge_simulations(force_field, md_name, sim_dirs):
  if not sim_dirs:
    return
  md_module = get_md_module(force_field)
  md_module.merge_simulations(md_name, sim_dirs)
  # calculate time spent in simulations
  fnames = ['%s/%s.time' % (pulse, md_name) 
            for pulse in sim_dirs
            if os.path.isfile(pulse)]
  vals = [open(f).read().split()[0] for f in fnames]
  time = sum([float(val) for val in vals])
  open(md_name+'.pulse.time', 'w').write(util.elapsed_time_str(time))

    
def pulse(
    force_field, in_name, md_name, n_step, pulse_fn, 
    n_step_per_pulse=100, restraint_pdb=""):
  """
  Takes as argument, a first order function:
    def pulse_fn(soup):
  that updates the velocities at the beginnig of every pulse
  """

  top, crds, vels = get_restart_files(in_name)

  parms = fetch_parms(
      force_field, top, crds, restraint_pdb, 'constant_energy', md_name)
  parms['input_md_name'] = in_name
  parms['input_vels'] = vels
  parms['n_step_dynamics'] = n_step
  parms['n_step_per_snapshot'] = n_step_per_pulse // 2
  parms['n_step_per_pulse'] = n_step_per_pulse

  config = md_name + ".config"
  if util.is_same_dict_in_file(parms, config):
    print "simulation already run."
    return

  timer = util.Timer()

  save_dir = os.getcwd()

  n_pulse = parms['n_step_dynamics'] / n_step_per_pulse
  n_step_list = [n_step_per_pulse for i in range(n_pulse)]
  n_excess_step = parms['n_step_dynamics'] % n_step_per_pulse
  if n_excess_step > 0:
    n_pulse += 1
    n_step_list.append(n_excess_step)
  
  pulse_parms = copy.deepcopy(parms)
  pulse_parms['topology'] = os.path.abspath(parms['topology'])
  if 'thermometer_type' in pulse_parms:
    del pulse_parms['thermometer_type']
  in_md_name = parms['input_md_name']
  pulse_parms['input_md_name'] = os.path.abspath(in_md_name)
  pulse_in_top, pulse_in_crds, pulse_in_vels = \
    get_restart_files(pulse_parms['input_md_name'])

  pulses = ["pulse%d" % i for i in range(n_pulse)]
  for pulse, n_step in zip(pulses, n_step_list):
    os.chdir(save_dir)

    print "Pulse: %s/%d" % (pulse, n_pulse)
    util.goto_dir(pulse)

    pulse_parms['n_step_dynamics'] = n_step

    soup = soup_from_restart_files(
        pulse_in_top, pulse_in_crds, pulse_in_vels)

    pulse_fn(soup)

    crds, vels = write_soup_to_crds_and_vels(
        force_field, soup, md_name + '.pulse.in')
    pulse_parms['input_crds'] = crds
    pulse_parms['input_vels'] = vels

    run_parms(pulse_parms)

    pulse_parms['input_md_name'] = os.path.abspath(md_name)

    pulse_in_top, pulse_in_crds, pulse_in_vels = \
      get_restart_files(md_name)

  os.chdir(save_dir)

  merge_simulations(force_field, md_name, pulses)
  util.clean_fname(*pulses)

  open(md_name+'.time', 'w').write(timer.str()+'\n')

  util.write_dict(config, parms)



    
    
