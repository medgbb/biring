"""
Author: Reza Salari
Name: biring.py  (Bond In Ring)
Version: 0.1
License: GPL v3

"""
from __future__ import (unicode_literals, print_function)
import os
import re
import argparse
from collections import defaultdict
import numpy as np

# =========================================================
RING_DEFINITIONS = """
# format:
# resname   atom names

PHE         CG  CZ  CD1 CD2 CE1 CE2
TYR         CG  CZ  CD1 CD2 CE1 CE2
HIS         CG  CD1 CD2 NE1 NE2
TRP         CG  CD1 CD2 NE1 CE2
TRP         CD2 CE2 CE3 CZ2 CZ3 CH2
CHL1        C1  C2  C3  C4  C5  C10
CHL1        C5  C6  C7  C8  C9  C10
CHL1        C8  C9  C11 C12 C13 C14
CHL1        C13 C14 C15 C16 C17

"""
# =========================================================


# =========================================================
# Helpers: ringdef parser, bond finder
# =========================================================
coord_dist = lambda crd1, crd2: np.dot(crd1 - crd2, crd1 - crd2) ** 0.5


def parse_ring_definitions(ringdefs):
    """ Parse the RING_DEFINITIONS.

    Returns a dictionary of the following format:
        parsed_rings: {"RESNAME":[ (RING1_ATOMS), (RING2_ATOMS),... ]}
        nres: int, number of residues
        nrings: int, number of rings
    """

    lines = ringdefs.split('\n')

    parsed_rings = defaultdict(list)
    nrings = 0

    for line in lines:
        line = line.strip()
        if line == '' or line.startswith('#'):
            continue
        fields = line.split()
        parsed_rings[fields[0]].append(tuple(fields[1:]))
        nrings += 1

    nres = len(list(parsed_rings.keys()))
    return parsed_rings, nres, nrings


def guess_intra_res_bonds(mol, target_resnames=[], bond_cutoff=3.0):
    """ Find the atoms that are bonded.

    This function uses an approximate method (geometric
        distance) to define bonds.

    Arguments:
        mol : a Molecule() instance
        target_resnames: list, ['POPC', 'DOPC', ...]
        bond_cutoff: float, cutoff for defining a bond

    Returns:
        resname_bonds_map = {'POPC':[ ('C1','C2'), ('C2','C3'), ... ], ...}
    """

    resname_bonds_map = {}

    for res in mol.residues:
        if res.name in target_resnames:
            if res.name not in list(resname_bonds_map.keys()):
                bonds = []
                N = res.natoms()
                for i in range(0, N - 1):
                    for j in range(i + 1, N):
                        crd1 = np.array(res.atoms[i].coord)
                        crd2 = np.array(res.atoms[j].coord)
                        dist = coord_dist(crd1, crd2)
                        if dist <= bond_cutoff:
                            bonds.append((res.atoms[i].name,
                                          res.atoms[j].name))
                resname_bonds_map[res.name] = bonds

    return resname_bonds_map


# =========================================================
# Classes
# =========================================================


class Atom:

    def __init__(self, name, number, index, coord):
        self.name = name
        self.number = number
        self.index = index
        self.coord = coord

        self.residue = None


class Residue:

    def __init__(self, name, number, index):
        self.name = name
        self.number = number
        self.index = index

        self.atoms = []
        self.molecule = None

    def add_atom(self, atname, atnumb, atindex, atcoord):
        a = Atom(atname, atnumb, atindex, atcoord)
        a.residue = self
        self.atoms.append(a)

        return a

    def get_atoms_basedon_names(self, *atnames):
        result = []
        for a in self.atoms:
            if a.name in atnames:
                result.append(a)

        return result

    def get_atom_coords_basedon_names(self, *atnames):
        result = []
        for a in self.atoms:
            if a.name in atnames:
                result.append(a.coord)

        return result

    def natoms(self):
        return len(self.atoms)

    def get_com(self):
        crds = [a.coord for a in self.atoms]
        return np.array(crds).mean(axis=0)

    def __repr__(self):
        return "%4s %4d [%4d]" % (self.name, self.number, self.index)


class Molecule:

    def __init__(self):
        self.residues = []

    def add_residue(self, resname, resnumber, resindex):
        r = Residue(resname, resnumber, resindex)
        r.molecule = self
        self.residues.append(r)

        return r

    def get_residues_by_names(self, *resnames):
        result = []
        for r in self.residues:
            if r.name in resnames:
                result.append(r)

        return result

    def nres(self):
        return len(self.residues)

    def natoms(self):
        return sum(r.natoms() for r in self.residues)

    def __repr__(self):
        return '(i) molecule with %d residues and %d atoms.' % (
            self.nres(), self.natoms())

# =========================================================
# PDB reader
# =========================================================


def read_pdb(fname, ignored_res=[], ignore_H=False):
    """ Read a PDB file.

    Arguments:
        fname: string, name of pdb file
        ignored_res: list, list of residue names to be ignored
            like ['TIP3',...]
        ignore_H: bool, whether to include the hydrogens in the molecule

    Returns:
        mol: a Molecule instance.

    """

    tmp_resname = None
    tmp_resnumb = None

    at_index = -1
    res_index = -1
    H_regex = re.compile('^[1-9]?[1-9]?H\w*')

    res = None

    mol = Molecule()

    with open(fname) as f:
        for line in f:
            line = line.strip()

            if (not line) or (not line.startswith(('ATOM', 'HETATM'))):
                continue

            atname = line[12:16].strip()
            resname = line[17:21].strip()
            resnumb = int(line[22:26].strip())

            x = line[30:38].strip()
            y = line[38:46].strip()
            z = line[46:54].strip()

            coord = tuple(map(float, [x, y, z]))

            at_index += 1

            if ignore_H:
                if atname[0] == 'H' or H_regex.match(atname):
                    continue

            if tmp_resname != resname or tmp_resnumb != resnumb:
                res_index += 1
                if resname in ignored_res:
                    continue
                res = mol.add_residue(resname, resnumb, res_index)
                tmp_resname = resname
                tmp_resnumb = resnumb

            res.add_atom(atname, None, at_index, coord)

    return mol

# =========================================================
# Function for finding the crossings
# =========================================================


def find_bonds_in_rings_alg1(mol, ring_defs, intruding_res, resname_bonds_map,
                             proximity_cutoff):
    """ Tries to find bonds inside rings using simple geometric distance.

    Arguments:
        mol: a Molecule() instance
        ring_defs: dict, output from parse_ring_definitions()
        intruding_res: list, list of residue names that are checked
            like ['POPC',...]
        resname_bonds_map: dict, output from guess_intra_res_bonds()
        proximity_cutoff: float, the cutoff between the center of the ring and
            the center of a bond for being reported.

    Returns:
        result: dict, {'res1-res2':[dist1, dist2,...], ....}

    """

    # find the index of residues that have rings
    ring_resnames = list(ring_defs.keys())
    ring_ind = []
    for i, res in enumerate(mol.residues):
        if res.name in ring_resnames:
            ring_ind.append((i, res.name))

    # find the index of residues that can pass through rings
    intrud_ind = []
    for j, res in enumerate(mol.residues):
        if res.name in intruding_res:
            intrud_ind.append((j, res.name))

    # result will be of format { 'CHL1 11 - POPC 34': [1.2, 1.4], ...}
    result = defaultdict(list)

    # result_res_objects will be of format [ (res1, res2), ...]
    result_res_objects = []

    for i, iresname in ring_ind:
        ires = mol.residues[i]
        for ring_atoms in ring_defs[iresname]:
            ring_crds = ires.get_atom_coords_basedon_names(*ring_atoms)
            ring_crds_mean = np.array(ring_crds).mean(axis=0)

            for j, jresname in intrud_ind:
                jres = mol.residues[j]
                # quick check to see if the ring and the intruder are near
                # each other
                jcom = jres.get_com()
                if coord_dist(ring_crds_mean, jcom) >= 20:
                    continue

                # now check the distance between the center of each bond to
                # the center of the ring
                for b in resname_bonds_map[jresname]:
                    bond_crds = jres.get_atom_coords_basedon_names(*b)
                    bond_crds_mean = np.array(bond_crds).mean(axis=0)
                    dist = coord_dist(ring_crds_mean, bond_crds_mean)
                    if dist <= 2.3:
                        key = '%s - %s' % (ires, jres)
                        result[key].append(dist)
                        if not (ires, jres) in result_res_objects:
                            result_res_objects.append((ires, jres))

    return result, result_res_objects

# =========================================================
# Make VMD script
# =========================================================


def make_vmd_script(fname, result_res_objects, vmd_script_name):
    """ Create a vmdscript for visualizing the result

    Arguments:
        fname: string, the name of PDB file
        result_res_objects: list, [(res1, res2), ...]

    """

    VMD_TEMPLATE = """
mol load pdb {FNAME}
mol delrep 0 top

{REP_COMMANDS}
    """

    script = VMD_TEMPLATE.replace('{FNAME}', fname)
    rep_commands = []

    rep_number = 0
    for ires, jres in result_res_objects:

        # sel = '(resname %s and resid %s) or (resname %s and resid %s)' % (
        #     ires.name, ires.number, jres.name, jres.number)
        sel = '(same residue as index %d) or (same residue as index %d)' % (
            ires.atoms[0].index, jres.atoms[0].index)
        cm = 'mol addrep top\n'
        cm += 'mol modselect %d top "noh and (%s)"\n' % (rep_number, sel)
        cm += 'mol showrep top %d 0\n' % (rep_number)

        rep_commands.append(cm)
        rep_number += 1

    script = script.replace('{REP_COMMANDS}', '\n'.join(rep_commands))

    if vmd_script_name is None:
        vmd_script_name = fname + '_vmd'
    with open(vmd_script_name, 'w') as f:
        f.writelines(script)

    return vmd_script_name

# =========================================================
# Main
# =========================================================


def main():
    # command line options
    p = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    p.add_argument('pdb', help='pdb file')

    p.add_argument('--cut', dest='bond_cutoff', default=1.6,
                   help='cutoff for defining bonds')

    p.add_argument('--prox', dest='proximity_cutoff', default=2.3,
                   help='proximity cutoff for finding bonds in rings')

    p.add_argument('-H', dest='include_H', default=False, action='store_true',
                   help='include hydrogens for checking bonds')

    p.add_argument('-I', dest='ignored_res', nargs='?', const=[],
                   default=['TIP3', 'TIP3P', 'SOD', 'CLA', 'HOH'],
                   help='ignore these residues to read the file faster')

    p.add_argument('-A', dest='intruding_res', nargs='?', const=[],
                   default=['POPC', 'DOPC'],
                   help='residues that can pass through a ring')

    p.add_argument('--vmd', dest='vmd_script_name', default=None,
                   help='the name of vmd script.')

    # args:
    # Namespace(bond_cutoff=2.8, include_H=False, pdb='step5_assembly.pdb')
    args = p.parse_args()

    # ---------------------------------
    # inform about flags
    print('(i) include hydrogens:', args.include_H)
    print('(i) cutoff for defining bonds: %4.2f A' % args.bond_cutoff)
    print('(i) proximity cutoff for reporting bonds in rings: %4.2f A' %
          args.proximity_cutoff)
    print('(i) ignoring these residues: %s' % args.ignored_res)
    print('(i) checking inruding residues: %s' % args.intruding_res)

    # ---------------------------------
    # check for the file
    if not os.path.exists(args.pdb):
        print('(e) the file "%s" doesn\'t exist' % args.pdb)
        return
    # read the molecule
    mol = read_pdb(args.pdb, args.ignored_res, not args.include_H)
    print(mol)

    # ---------------------------------
    # guess bonds
    resname_bonds_map = guess_intra_res_bonds(mol, args.intruding_res,
                                              args.bond_cutoff)
    print('(i) number of found bonds based on %4.2f A cutoff: ' %
          args.bond_cutoff)
    for k in list(resname_bonds_map.keys()):
        print('\t %6s : %2d bonds' % (k, len(resname_bonds_map[k])))

    # ---------------------------------
    # parse ring definitons
    ring_defs, nres, nrings = parse_ring_definitions(RING_DEFINITIONS)
    print('(i) parsed %d rings belonging to %d residues.' % (nrings, nres))

    # ---------------------------------
    # find suspicious residues using algorithm 1
    print('(i) finding bonds in rings...')
    result, result_res_objects = \
        find_bonds_in_rings_alg1(mol,
                                 ring_defs,
                                 args.intruding_res,
                                 resname_bonds_map,
                                 args.proximity_cutoff)

    # show the result - ascending based on the distance
    print('(i) these residues potentially have bond-in-ring situation:\n')
    result_sum = {k: min(result[k]) for k in list(result.keys())}
    result_keys = sorted(result_sum, key=result_sum.get)
    for k in result_keys:
        print('\t %s    -> %4.2f A' % (k, result_sum[k]))
    print(' ')

    # ---------------------------------
    # make vmd script
    outname = make_vmd_script(args.pdb, result_res_objects,
                              args.vmd_script_name)
    print('(i) the vmd script is created - to visulize the reported residues:')
    print('\n         vmd -e %s\n' % outname)
    print('    note that the representations are hidden by default.')
    print('    go to Graphics->Representions to enable them.\n')

    # ---------------------------------
    print('(i) done.')


if __name__ == '__main__':
    main()
