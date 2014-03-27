
# Biring

Reza Salari - [Brannigan Lab](http://branniganlab.org)


## Introduction

When setting up a simulation box (especially when lipids are involved), it is
possible that some carbon chains end up passing through rings by mistake (e.g.
tails of POPC goes through a benzene ring of PHE). This script tries to find
such cases.

Currently, one simple algorithm is implemented - if the distance between the
center of a bond from the center of a ring is less than the `PROXIMITY_CUTOFF`,
it will be reported. The script also creates a vmd script that help visualizing
the reported residues.

To make the script flexible and simple, it only needs a PDB file and doesn't
rely on the topology files (like PSF or TOP). On the flip side, it uses the
cutoff method to find the atoms that form bonds. The current `BOND_CUTOFF` is
1.6 A (a bit longer than the length of C-C bond).

## Requirements

- Python 2.7 or higher
- Numpy

## Defining new rings
You can add more ring definitions to the `RING_DEFINITIONS` in the beginning of
the script.

## How to use it

Use `-h` for all options. Typically you can run:

    python biring.py pdbfname.pdb

And it will output the residues that it thinks might have bond-in-ring issue.

Sample output:

	$ python3.4 biring.py sample/pdb1.pdb

	(i) include hydrogens: False
	(i) cutoff for defining bonds: 1.60 A
	(i) proximity cutoff for reporting bonds in rings: 2.30 A
	(i) ignoring these residues: ['TIP3', 'TIP3P', 'SOD', 'CLA', 'HOH']
	(i) checking inruding residues: ['POPC', 'DOPC']
	(i) molecule with 1919 residues and 23471 atoms.
	(i) number of found bonds based on 1.60 A cutoff:
	           POPC : 49 bonds
	(i) parsed 9 rings belonging to 5 residues.
	(i) finding bonds in rings...
	(i) these residues potentially have bond-in-ring situation:

          PHE  328 [1007] - POPC   40 [1738]    -> 1.14 A
         CHL1   17 [1715] - POPC   35 [1733]    -> 1.30 A
         CHL1    1 [1699] - POPC   62 [1760]    -> 1.36 A
         CHL1   29 [1727] - POPC   39 [1737]    -> 1.36 A
          TRP  308 [ 987] - POPC  148 [1846]    -> 1.64 A
         CHL1  122 [1820] - POPC  137 [1835]    -> 1.80 A
          PHE  327 [ 326] - POPC   45 [1743]    -> 1.85 A
         CHL1   24 [1722] - POPC   91 [1789]    -> 1.92 A
         CHL1  108 [1806] - POPC  216 [1914]    -> 2.07 A
         CHL1   26 [1724] - POPC   79 [1777]    -> 2.17 A
         CHL1   28 [1726] - POPC   88 [1786]    -> 2.25 A
         CHL1  123 [1821] - POPC  149 [1847]    -> 2.30 A

	(i) the vmd script is created - to visulize the reported residues:

	         vmd -e sample/pdb1.pdb_vmd

	    note that the representations are hidden by default.
	    go to Graphics->Representions to enable them.

	(i) done.

Note that the results need visual inspection, as the current algorithm is simple
and a bit paranoid. In the case of the sample above, only two sets of residues
(PHE 328 - POPC 40) and (CHL1 122 - POPC 137) indeed have bonds crossing the
rings.

## To run tests

	python test_biring.py

	# or if you have coverage installed:
	python -m coverage run test_biring.py
	python -m coverage report
	python -m coverage html

## Questions

Please use the [issues](https://github.com/resal81/biring/issues) page.
