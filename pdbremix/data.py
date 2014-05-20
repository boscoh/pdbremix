# encoding: utf-8

__doc__ = """ 

Stores data for PDBREMIX.

- directory of explici data directory
- mappings of residue names to chars
- names and locations of binaries used
- backbone atom types
- solvent residue types
- element radii
- chi dihedral angle topologies
"""


import os
import util


module_dir = os.path.dirname(__file__)
data_dir = os.path.join(module_dir, 'data')

binaries_fname = os.path.join(data_dir, 'binaries.json')
binaries = util.read_dict(binaries_fname)

def binary(bin, arg_str='', out_name=None, in_fname=None):
  """
  Runs an external binary, handles arguments, writes out
  equivalent .sh file, log file, and can pipe in in_fname.
  """
  if bin in binaries and binaries[bin]:
    bin = binaries[bin]
  else:
    util.check_program(bin)
  if arg_str:
    util.run_with_output_file(
        '%s %s' % (bin, arg_str), out_name, in_fname)
  return '"%s"' % bin


def invert_dict(d):
  """
  Returns dictionaries with swapped key-values.
  """
  return dict((v, k) for k,v in d.items())


res_name_to_char = {
    "ALA":"A", "CYS":"C", "ASP":"D",
    "GLU":"E", "PHE":"F", "GLY":"G",
    "HIS":"H", "ILE":"I", "LYS":"K",
    "LEU":"L", "MET":"M", "ASN":"N",
    "PRO":"P", "GLN":"Q", "ARG":"R",
    "SER":"S", "THR":"T", "VAL":"V",
    "TRP":"W", "TYR":"Y", "ACE":">",
    "NME":"<",
}

res_char_to_name = invert_dict(res_name_to_char)


# recognized atom types for protein backbone 
backbone_atoms = [
    "OXT", # C-terminal carboxyl group
    "H1", "H2", "H3",  # N-terminal charged group
    "C", "O", # peptide-bond carbonyl group
    "H", "HN", "N",  # peptide-bond amide group
    "CA", "HA" # main-chain C-alpha and alkyl group
    ]
    
solvent_res_types = [
    'HOH', 'WAT', 'TIP', 'SOL',
    'CLA', 'SOD', 'NA', 'CL', 
    'NA+', 'CL-', 'Na', 'Cl',
    'Na+', 'Cl-']

radii = { 
 'H':  1.20,
 'N':  1.55,
 'NA': 2.27,
 'CU': 1.40,
 'CL': 1.75,
 'C':  1.70,
 'O':  1.52,
 'I':  1.98,
 'P':  1.80,
 'B':  1.85,
 'BR': 1.85,
 'S':  1.80,
 'SE': 1.90,
 'F':  1.47,
 'FE': 1.80,
 'K':  2.75,
 'MN': 1.73,
 'MG': 1.73,
 'ZN': 1.39,
 'HG': 1.80,
 'XE': 1.80,
 'AU': 1.80,
 'LI': 1.80,
 '.':  1.80
}

two_char_elements = [e for e in radii.keys() if len(e) == 2]


def strip_numbers(s):
  result = ""
  for c in s:
    if not c.isdigit() and c != " ":
      result += c
  return result


def guess_element(res_type, atom_type):
  """
  Returns the element type using a dirty heuristic guess.
  """
  atom_type = strip_numbers(atom_type)
  if res_type in res_name_to_char:
    return atom_type[0]
  if len(atom_type) == 2 and atom_type in two_char_elements:
    return atom_type
  return atom_type[0]  


chi_topology = {
  'ARG': [ ['N', 'CA', 'CB', 'CG'],
           ['CA', 'CB', 'CG', 'CD'],
           ['CB', 'CG', 'CD', 'NE'],
           ['CG', 'CD', 'NE', 'CZ']],
  'ASN': [['N', 'CA', 'CB', 'CG'], ['CA', 'CB', 'CG', 'OD1']],
  'ASP': [['N', 'CA', 'CB', 'CG'], ['CA', 'CB', 'CG', 'OD1']],
  'CYS': [['N', 'CA', 'CB', 'SG']],
  'GLN': [ ['N', 'CA', 'CB', 'CG'],
           ['CA', 'CB', 'CG', 'CD'],
           ['CB', 'CG', 'CD', 'OE1']],
  'GLU': [ ['N', 'CA', 'CB', 'CG'],
           ['CA', 'CB', 'CG', 'CD'],
           ['CB', 'CG', 'CD', 'OE1']],
  'HIS': [['N', 'CA', 'CB', 'CG'], ['CA', 'CB', 'CG', 'ND1']],
  'ILE': [['N', 'CA', 'CB', 'CG1'], ['CA', 'CB', 'CG1', 'CD1']],
  'LEU': [['N', 'CA', 'CB', 'CG'], ['CA', 'CB', 'CG', 'CD1']],
  'LYN': [ ['N', 'CA', 'CB', 'CG'],
           ['CA', 'CB', 'CG', 'CD'],
           ['CB', 'CG', 'CD', 'CE'],
           ['CG', 'CD', 'CE', 'NZ']],
  'LYP': [ ['N', 'CA', 'CB', 'CG'],
           ['CA', 'CB', 'CG', 'CD'],
           ['CB', 'CG', 'CD', 'CE'],
           ['CG', 'CD', 'CE', 'NZ']],
  'LYS': [ ['N', 'CA', 'CB', 'CG'],
           ['CA', 'CB', 'CG', 'CD'],
           ['CB', 'CG', 'CD', 'CE'],
           ['CG', 'CD', 'CE', 'NZ']],
  'MET': [ ['N', 'CA', 'CB', 'CG'],
           ['CA', 'CB', 'CG', 'SD'],
           ['CB', 'CG', 'SD', 'CE']],
  'PHD': [['N', 'CA', 'CB', 'CG'], ['CA', 'CB', 'CG', 'OD1']],
  'PHE': [['N', 'CA', 'CB', 'CG'], ['CA', 'CB', 'CG', 'CD1']],
  'PRO': [ ['N', 'CA', 'CB', 'CG'],
           ['CA', 'CB', 'CG', 'CD'],
           ['CB', 'CG', 'CD', 'N'],
           ['CG', 'CD', 'N', 'CA']],
  'SER': [['N', 'CA', 'CB', 'OG']],
  'THR': [['N', 'CA', 'CB', 'OG1']],
  'TRP': [['N', 'CA', 'CB', 'CG'], ['CA', 'CB', 'CG', 'CD1']],
  'TYR': [['N', 'CA', 'CB', 'CG'], ['CA', 'CB', 'CG', 'CD1']],
  'VAL': [['N', 'CA', 'CB', 'CG1']]}


def get_res_chi_topology(res_type):
  """
  Returns the chi topology for a given residue, which is a list of
  atoms that are affected if one rotates the chi0, chi1... 
  dihedral angle.
  """
  # Some common residue renamings in AMBER and GROMACS
  if res_type in ["HID", "HIE", "HIP", "HSE"] or "HIS" in res_type:
    res_type = "HIS"
  if res_type in ["LYP"]:
    res_type = "LYS"
  if res_type in ["CYM", "CYX", "CYN"]:
    res_type = "CYS"
  if res_type not in chi_topology:
    return []
  else:
    return chi_topology[res_type]



