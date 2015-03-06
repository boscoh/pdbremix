import array
import math
import v3


class SpaceHash(object):

    def __init__(self, vertices, div=5.3, padding=0.05):
        self.vertices = vertices
        self.div = div
        self.inv_div = 1.0 / self.div
        self.padding = padding

        zero3 = lambda: [0.0] * 3
        self.minima = zero3()
        self.maxima = zero3()
        self.spans = zero3()
        self.sizes = zero3()

        for i in range(3):
            self.minima[i] = min([v[i] for v in self.vertices])
            self.maxima[i] = max([v[i] for v in self.vertices])
            self.minima[i] -= self.padding
            self.maxima[i] += self.padding
            self.spans[i] = self.maxima[i] - self.minima[i]
            self.sizes[i] = int(math.ceil(self.spans[i] * self.inv_div))

        self.size1_size2 = self.sizes[1] * self.sizes[2]
        self.size2 = self.sizes[2]

        self.cells = {}
        self.spaces = []
        for i_vertex, vertex in enumerate(self.vertices):
            space = self.vertex_to_space(vertex)
            self.spaces.append(space)
            space_hash = self.space_to_hash(space)
            cell = self.cells.setdefault(space_hash, array.array('L'))
            cell.append(i_vertex)

    def vertex_to_space(self, v):
        return [int((v[i] - self.minima[i]) * self.inv_div) for i in range(3)]

    def space_to_hash(self, s):
        return s[0] * self.size1_size2 + s[1] * self.size2 + s[2]

    def neighbourhood(self, space):
        def neighbourhood_in_dim(space, i_dim):
            i = max(0, space[i_dim] - 1)
            j = min(self.sizes[i_dim], space[i_dim] + 2)
            return range(i, j)
        for s0 in neighbourhood_in_dim(space, 0):
            for s1 in neighbourhood_in_dim(space, 1):
                for s2 in neighbourhood_in_dim(space, 2):
                    yield [s0, s1, s2]

    def close_pairs(self):
        n_vertex = len(self.vertices)
        for i_vertex0 in range(n_vertex):
            space0 = self.spaces[i_vertex0]
            for space1 in self.neighbourhood(space0):
                hash1 = self.space_to_hash(space1)
                for i_vertex1 in self.cells.get(hash1, []):
                    if i_vertex0 < i_vertex1:
                        yield i_vertex0, i_vertex1


def find_bb_hbonds(residues):
    vertices = []
    atoms = []
    for i_residue, residue in enumerate(residues):
        residue.i = i_residue
        for atom in residue.atoms():
            atom.residue = residue
        if residue.has_atom('O'):
            atom = residue.atom('O')
            atoms.append(atom)
            vertices.append(atom.pos)
        if residue.has_atom('N'):
            atom = residue.atom('N')
            atoms.append(atom)
            vertices.append(atom.pos)
    d = 3.5
    for i, j in SpaceHash(vertices).close_pairs():
        if atoms[i].element == 'O' and atoms[j].element == 'N':
            o = atoms[i]
            n = atoms[j]
        elif atoms[i].element == 'N' and atoms[j].element == 'O':
            n = atoms[i]
            o = atoms[j]
        else:
            continue
        if abs(o.residue.i - n.residue.i) < 3:
            continue
        if v3.distance(o.pos, n.pos) < d:
            o.residue.co_partners.append(n.residue.i)
            n.residue.nh_partners.append(o.residue.i)


def find_ss_by_bb_hbonds(soup):
    residues = soup.residues()
    n_res = len(residues)

    for res in residues:
        res.co_partners = []
        res.nh_partners = []
        res.beta_contacts = []
        res.alpha_contacts = []
        res.ss = "C"

    find_bb_hbonds(residues)

    def is_conh(i_res, j_res):
        if not (0 <= i_res < n_res):
            return False
        if not (0 <= j_res < n_res):
            return False
        return j_res in residues[i_res].co_partners

    def unique_append(a_list, item):
        if item not in a_list:
            a_list.append(item)
            a_list.sort()

    def make_alpha_contacts(i_res, j_res):
        unique_append(residues[i_res].alpha_contacts, j_res)
        unique_append(residues[j_res].alpha_contacts, i_res)

    for i_res in range(n_res):
        # alpha-helix
        if is_conh(i_res - 1, i_res + 3) and is_conh(i_res, i_res + 4):
            for j_res in range(i_res, i_res + 4):
                residues[j_res].ss = 'H'
            make_alpha_contacts(i_res + 3, i_res)
            make_alpha_contacts(i_res + 4, i_res)

        # 3-10 helix
        if is_conh(i_res - 1, i_res + 2) and is_conh(i_res, i_res + 3):
            for j_res in range(i_res, i_res + 3):
                residues[j_res].ss = 'H'
            make_alpha_contacts(i_res + 3, i_res)

    def make_beta_contacts(i_res, j_res):
        unique_append(residues[i_res].beta_contacts, j_res)
        unique_append(residues[j_res].beta_contacts, i_res)

    for i_res in range(n_res):
        for j_res in range(n_res):
            if abs(i_res - j_res) < 3:
                continue

            # anti-parallel beta-sheet h-bonded pair
            if is_conh(i_res, j_res) and is_conh(j_res, i_res):
                make_beta_contacts(i_res, j_res)
                residues[i_res].ss = "E"
                residues[j_res].ss = "E"

            # anti-parallel beta-sheet non-h-bonded pair
            if is_conh(i_res - 1, j_res + 1) and is_conh(i_res + 1, j_res - 1):
                residues[i_res].ss = "E"
                residues[j_res].ss = "E"
                make_beta_contacts(i_res, j_res)

            # parallel beta sheet pairs
            if is_conh(i_res, j_res - 1) and is_conh(j_res - 1, i_res):
                make_beta_contacts(i_res, j_res)
                residues[i_res].ss = "E"
                residues[j_res].ss = "E"
    
    # TODO: check for anti-parallel cross-strand i->j-2 contacts

def is_ss_contact(soup, i_res, j_res):
    residue = soup.residues()[i_res]
    return j_res in residue.beta_contacts or \
           j_res in residue.alpha_contacts

