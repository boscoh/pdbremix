import v3
import copy
import string

import data

class Atom:
  def __init__(
      self, pos=None,
      atom_type="", res_num=None):
    self.is_hetatm = False
    self.pos = v3.vector() if pos is None else pos
    self.vel = v3.vector()
    self.mass = 0.0
    self.type = ""
    self.element = ""
    self.chain_id = " "
    self.res_type = ""
    self.res_num = ""
    self.res_insert = ""
    self.bfactor = 0.0
    self.occupancy = 0.0
    self.num = 0
    self.alt_conform = " "

  def copy(self):
    return copy.deepcopy(self)

  def type_str(self):
    atom_type = self.type
    if len(atom_type) == 1:
      atom_type = " %s  " % atom_type
    elif len(atom_type) == 2:
      atom_type = " %s " % atom_type
    elif len(atom_type) == 3:
      if atom_type[0].isdigit():
        atom_type = "%s " % atom_type
      else:
        atom_type = " %s" % atom_type
    return atom_type    

  def pdb_str(self):
    if self.is_hetatm:
      field = "HETATM"
    else:
      field = "ATOM  "
    x, y, z = self.pos
    s = "%6s%5s %4s %-4s%1s%4s%1s   %8.3f%8.3f%8.3f%6.2f%6.2f" \
            % (field, 
               str(self.num)[-5:], 
               self.type_str(),
               self.res_type, 
               self.chain_id,
               str(self.res_num)[-4:], 
               self.res_insert,
               x, y, z,
               self.occupancy, 
               self.bfactor)
    return s
               
  def __str__(self):
    x, y, z = self.pos
    return "%s%s-%s (% .1f % .1f % .1f)" \
            %  (self.res_type, self.res_num, 
                self.type, x, y, z)

  def transform(self, matrix):
    new_pos = v3.transform(matrix, self.pos)
    v3.set_vector(self.pos, new_pos)


def AtomFromPdbLine(line):
  """Returns an Atom object from an atom line in a pdb file."""
  atom = Atom()
  if line.startswith('HETATM'):
    atom.is_hetatm = True
  else:
    atom.is_hetatm = False
  atom.num = int(line[6:11])
  atom.type = line[12:16].strip(" ")
  element = ''
  for c in line[12:15]:
    if not c.isdigit() and c != " ":
      element += c
  if element[:2] in two_char_elements:
    atom.element = element[:2]
  else:
    atom.element = element[0]
  atom.alt_conform = line[16]
  atom.res_type = line[17:21].strip()
  atom.chain_id = line[21]
  atom.res_num = int(line[22:26])
  atom.res_insert = line[26]
  if atom.res_insert == " ":
    atom.res_insert = ""
  x = float(line[30:38])
  y = float(line[38:46])
  z = float(line[46:54])
  v3.set_vector(atom.pos, x, y, z)
  try:
    atom.occupancy = float(line[54:60])
  except:
    atom.occupancy = 100.0
  try:
    atom.bfactor = float(line[60:66])
  except:
    atom.bfactor = 0.0
  return atom
  
  
def cmp_atom(a1, a2):
  if a1.num < a2.num:
    return -1
  else:
    return 0


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


def guess_element(res_type, atom_type):
  if res_type in data.res_name_to_char:
    return atom_type[0]
  element = ""
  for c in atom_type:
    if not c.isdigit() and c != " ":
      element += c
  if len(element) == 2 and element in two_char_elements:
    return element
  return element[0]  


def add_radii(atoms):
  for atom in atoms:
    if atom.element in radii:
      atom.radius = radii[atom.element]
    else:
      atom.radius = radii['.']


def get_center(atoms):
  center = v3.vector()
  for atom in atoms:
    center += atom.pos
  result = v3.scale(center, 1.0/float(len(atoms)))
  return result


def get_width(atoms, center):
  max_diff = 0
  for atom in atoms:
    diff = v3.distance(atom.pos, center)
    if diff > max_diff:
      max_diff = diff
  return 2*max_diff


class AtomList:

  def __init__(self, pdb=""):
    self.id = ''
    self._atoms = []
    if pdb:
      self.read_pdb(pdb)

  def copy(self):
    return copy.deepcopy(self)

  def n_atom(self):
    return len(self._atoms)

  def atoms(self):
    return self._atoms

  def atom(self, i):
    return _atoms[i]
    
  def clear(self):
    for atom in self._atoms:
      del atom
    del self._atoms[:]

  def transform(self, matrix):
    for atom in self._atoms:
      atom.transform(matrix)

  def insert_atom(self, atom):
    self._atoms.append(atom)
    
  def erase_atom(self, atom_type):
    for atom in self._atoms:
      if atom.type == atom_type:
        self._atoms.remove(atom)
        del atom
        return

  def read_pdb(self, fname):
    self.clear()
    for line in open(fname, 'r').readlines():
      if line.startswith(("ATOM", "HETATM")):
        atom = AtomFromPdbLine(line);
        if len(self._atoms) == 1:
          self.id = atom.chain_id
        self.insert_atom(atom)
      if line.startswith("ENDMDL"):
        return

  def write_pdb(self, pdb):
    with open(pdb, 'w') as f:
      for atom in sorted(self._atoms, cmp=cmp_atom):
        f.write(atom.pdb_str() + '\n')

  def set_id(self, new_id):
    self.id = new_id
    for a in self.atoms():
      a.chain_id = new_id



class Residue:

  def __init__(self, in_type, in_chain_id, in_num, in_insert=''):
    self.type = in_type
    self.chain_id = in_chain_id
    self.num = in_num
    self.insert = in_insert
    self._atom_dict = {}
 
  def name(self):
    tag = ""
    if self.chain_id != " " and self.chain_id != "":
      tag += self.chain_id + ":"
    tag += str(self.num)
    if self.insert:
      tag += self.insert
    return tag  

  def __str__(self):
    atom_name_list = [a.type for a in self.atoms()]
    atom_name = " ".join(atom_name_list)
    return "%s-%s { %s }" % (self.type, self.num, atom_name)

  def copy(self):
    return copy.deepcopy(self)
  
  def n_atom(self):
    return len(self._atom_dict)
    
  def atom(self, atom_type):
    return self._atom_dict[atom_type]
    
  def has_atom(self, atom_type):
    return atom_type in self._atom_dict.keys()
    
  def change_atom_type(self, atom_type1, atom_type2):
    atom = self._atom_dict[atom_type1]
    atom.type = atom_type2
    del self._atom_dict[atom_type1]
    self._atom_dict[atom_type2] = atom

  def atoms(self):
    return self._atom_dict.values()
  
  def atom_name(self, atom_type):
    return self.type + self.num + ":" + atom_type

  def insert_atom(self, atom):
    self._atom_dict[atom.type] = atom
    atom.chain_id = self.chain_id
    atom.res_num = self.num
    atom.res_type = self.type
  
  def erase_atom(self, atom_type):
    del self._atom_dict[atom_type]
    
  def set_num(self, i, insert=""):
    self.num = i
    self.insert = insert
    for atom in self.atoms():
      atom.res_num = self.num
      atom.res_insert = insert
     
  def inc_num(self):
    self.set_num(self.num+1, self.insert)

  def dec_num(self):
    self.set_num(self.num-1, self.insert)
    
  def dec_insert(self):
    l = self.insert;
    if l == "A" or l == "a":
      self.insert = ''
    else:
      i = string.ascii_letters.find(l)
      self.insert = string.ascii_letters[i-1]

  def transform(self, matrix):
     for atom in self.atoms():
       atom.transform(matrix)

  def set_chain_id(self, chain_id):
    self.chain_id = chain_id
    for a in self.atoms():
      a.chain_id = chain_id

  def set_type(self, res_type):
    self.type = res_type
    for a in self.atoms():
      a.res_type = res_type


class Polymer(AtomList):

  def __init__(self, fname=""):
    AtomList.__init__(self)
    self._residues = []
    if fname:
      self.read_pdb(fname)

  def copy(self):
    return copy.deepcopy(self)

  def residue(self, i):
    return self._residues[i]
    
  def residues(self):
    return self._residues

  def insert_atom(self, i, atom):
    self._atoms.append(atom)
    self.residue(i).insert_atom(atom)
    
  def erase_atom(self, i, atom_type):
    atom = self.residue(i).atom(atom_type)
    self._atoms.remove(atom)
    self.residue(i).erase_atom(atom_type)
    del atom
    
  def clear(self):
    del self._residues[:]
    AtomList.clear(self)
    
  def n_residue(self):
    return len(self._residues)
    
  def insert_residue(self, i, res):
    is_insertion = False
    if i < self.n_residue()-1:
      save_res_num = self.residue(i).num
      if self.residue(i+1).num == save_res_num:
        is_insertion = True

    if self.n_residue() == 0:
      res.set_num(res.num, res.insert)
    elif i < self.n_residue():
      res.set_num(self.residue(i).num, self.residue(i).insert)
    else:
      res.set_num(self.residue(i-1).num, "")
      res.inc_num()

    self._residues.insert(i, res)
    for atom in res.atoms():
      self.insert_atom(i, atom)

    for j in range(i+1, self.n_residue()):
      self.residue(j).inc_num()

    if is_insertion:
      while self.residue(i+1).insert:
        for j in range(i+1, self.n_residue()):
          if self.residue(j).res_num == save_res_num:
            self.residue(k).dec_insert()
    
  def append_residue(self, res):
    self._residues.append(res)
    for atom in res.atoms():
      self.insert_atom(self.n_residue()-1, atom)

  def erase_residue(self, i):  
    save_res_num = self.residue(i).num

    for atom in self.residue(i).atoms():
      self._atoms.remove(atom)
      del atom
    self._residues.pop(i)  
    
    if i < self.n_residue():
      if self.residue(i).num == save_res_num:
        # erasing residue in an insertion
        for j in range(i, self.n_residue()):
          if self.residue(j).num == erase_res_num_int:
            self.residue(j).dec_insert()
      else:
        for j in range(i, self.n_residue()):
          self.residue(j).dec_num()
    
  def get_i_residue(self, tag):
    # clean up tag
    tag = tag.strip()
    if tag[0] == ":":
      tag = tag[1:]
    if not tag[0].isdigit() and tag[1].isdigit():
      tag = tag[0] + ":" + tag[1:]
    for i, residue in enumerate(self.residues()):
      if tag.lower() == residue.tag().lower():
        return i
    raise Exception("Can't find residue " + tag)
  
  def extract_polymer(self, i, j):
    extract = Polymer()
    for res in self.residues()[i:j]:
      extract.append_residue(res.copy())
    return extract
 
  def chain_ids(self):
    chain_id = [r.chain_id for r in self.residues()]
    return list(set(chain_id))

  def extract_chain(self, chain_id):
    extract = Polymer()
    for res in self.residues():
      if res.chain_id == chain_id:
        extract.append_residue(res.copy())
    return extract
 
  def insert_polymer(self, i, insert):
    for res in reversed(insert.residues()):
      self.insert_residue(i, res.copy())
    
  def load_residue_bfactors(self, res_bfactors):
    for i, r in enumerate(self.residues()):
      for atom in r.atoms():
        if i >= len(res_bfactors):
          return
        else:
          atom.bfactor = res_bfactors[i]

  def __str__(self):
    res_name_list = [str(res) for res in self._residues]
    return "\n".join(res_name_list)
 
  def read_pdb(self, fname):
    self.clear()
    res_num = -1
    res_insert = " "
    for line in open(fname, 'r').readlines():
      if line.startswith("ATOM") or line.startswith("HETATM"):
        atom = AtomFromPdbLine(line);
        if (res_num != atom.res_num) or (res_insert != atom.res_insert):
          residue = Residue(atom.res_type, atom.chain_id,
                            atom.res_num, atom.res_insert)
          self.append_residue(residue)
          res_num = atom.res_num
          res_insert = atom.res_insert
        self.insert_atom(-1, atom)
      if line.startswith("ENDMDL"):
        return

  def write_pdb(self, pdb):
    f = open(pdb, 'w')
    n_atom = 0
    for res in self.residues():
      res_atoms = res.atoms()
      res_atoms.sort(cmp_atom)
      for atom in res_atoms:
        f.write(atom.pdb_str() + '\n')
    f.close()
