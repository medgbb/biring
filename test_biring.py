
import unittest
import biring
import numpy as np

SAMPLE_RING_DEFINITIONS = """
# format:
# resname   atom names

PHE         CG  CZ  CD1 CD2 CE1 CE2
TYR         CG  CZ  CD1 CD2 CE1 CE2
CHL1        C1  C2  C3  C4  C5  C10
CHL1        C8  C9  C11 C12 C13 C14
CHL1        C13 C14 C15 C16 C17

"""


class TestBiringHelpers(unittest.TestCase):

    def test_coord_dist(self):
        c1 = np.array([0, 0, 0])
        c2 = np.array([0, 0, 5])
        self.assertEqual(biring.coord_dist(c1, c2), 5)

        c3 = np.array([4, 6, 7])
        result = np.sqrt(4. ** 2 + 6. ** 2 + 7. ** 2)
        self.assertEqual(biring.coord_dist(c1, c3), result)

    def test_parse_ring_definitions(self):
        parsed_defs, nres, nrings = \
            biring.parse_ring_definitions(SAMPLE_RING_DEFINITIONS)
        self.assertEqual(nres, 3)
        self.assertEqual(nrings, 5)
        self.assertEqual(
            parsed_defs['CHL1'][0], ('C1', 'C2', 'C3', 'C4', 'C5', 'C10'))

    def test_guess_intra_res_bonds(self):
        pass


class TestBiringMolClasses(unittest.TestCase):

    def test_residue(self):
        r = biring.Residue('SER', 1, 0)

        self.assertEqual(len(r.atoms), r.natoms())
        self.assertEqual(r.natoms(), 0)

        r.add_atom('N', None, 1, (0.1, 3.4, 9.0))
        self.assertEqual(r.natoms(), 1)
        a = r.atoms[0]
        self.assertEqual(a.name, 'N')
        self.assertEqual(a.coord, (0.1, 3.4, 9.0))
        self.assertEqual(a.residue, r)

        r.add_atom('CA', None, 2, (1.1, 2.2, -3.0))
        self.assertEqual(r.natoms(), 2)
        com = [round(x, 2) for x in r.get_com()]
        self.assertEqual(com, [0.6, 2.8, 3.0])


class TestBiringPDBReader(unittest.TestCase):

    def test_pdb_reader(self):
        mol = biring.read_pdb('sample/pdb1.pdb')

        self.assertEqual(mol.natoms(), 52803)
        self.assertEqual(mol.nres(), 1919)

        self.assertEqual(mol.residues[0].name, 'SER')
        self.assertEqual(mol.residues[-1].name, 'POPC')

        self.assertEqual(mol.residues[0].number, 1)
        self.assertEqual(mol.residues[-1].number, 220)

        self.assertEqual(mol.residues[0].atoms[0].name, 'N')
        self.assertEqual(mol.residues[-1].atoms[-1].name, 'H16Z')

        self.assertEqual(mol.residues[0].atoms[0].index, 0)
        self.assertEqual(mol.residues[-1].atoms[-1].index, 52802)


class TestBiringFindBondsInRingsAlg1(unittest.TestCase):
    pass

if __name__ == '__main__':
    unittest.main()
