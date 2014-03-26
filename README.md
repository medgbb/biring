
Reza Salari ([Brannigan Lab](http://branniganlab.org))


## Introduction

When setting up a simulation box (especially when lipids are involved), it is
possible that some carbon chains end up passing through rings by mistake (e.g.
tails of POPC goes through a benzene ring of PHE). This script tries to find
such cases.

Currently, one simple algorithm is implemented - if the distance between the
center of a bond from the center of a ring is less than the `PROXIMITY_CUTOFF`,
it will be reported. The script also creates a vmd script that help visualizing
the reported residues.

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


## To run tests

	python test_biring.py

	# or if you have coverage installed:
	python -m coverage run test_biring.py
	python -m coverage report
