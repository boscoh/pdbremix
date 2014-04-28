import os
import util

module_dir = os.path.dirname(__file__)
data_dir = os.path.join(module_dir, 'data')

def invert_dict(d):
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

binaries_fname = os.path.join(data_dir, 'binaries.json')
binaries = util.read_dict(binaries_fname)

def binary(bin, arg_str='', out_name=None, in_fname=None):
  if bin in binaries:
    bin = binaries[bin]
  util.check_program(bin)
  if arg_str:
    util.run_with_output_file('"%s" %s' % (bin, arg_str), out_name, in_fname)
  return '"%s"' % bin

# recognized atom types for protein backbone in AMBER
backbone_atoms = [
    "OXT", # C-terminal carboxyl group
    "H1", "H2", "H3",  # N-terminal charged group
    "C", "O", # peptide-bond carbonyl group
    "H", "HN", "N",  # peptide-bond amide group
    "CA", "HA" # main-chain C-alpha and alkyl group
    ]
    
solvent_res_types = [
    'WAT', 'TIP', 'HOH', 'SOL', 'CLA', 'SOD',
    'NA', 'CL', 'NA+', 'CL-']