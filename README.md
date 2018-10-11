# MolViewer
**MolViewer** is a JavaScript library for drawing Chemical Structures. It brings tools to easily convert MOL files and SMILE strings into 2D and 3D representations.

## Resources
* [API Reference](#api)
* [Examples](https://github.com/BoboRett/MolViewer/Examples)


## API

### Classes
#### Atom( *int* index, *array* position, *string* element[, *int* charge = 0] )
#### Bond( *int* index, *Atom* bondStart, *Atom* bondEnd[, *int* bondType = 1, *int* bondDirection = 0] )
#### Molecule( [*string* molFile] )

##### Properties
* *array[[Atom][Atom]]* Molecule.**atoms**

    List of molecule's atoms.

* *array[[Bond][Bond]]* Molecule.**bonds**

    List of molecule's bonds.

* *array* Molecule.**fGroups**

    List of molecule's functional groups. Requires Molecule.fGroupSearcher() to populate.

##### Methods

* Molecule.**parseMol**( *string* molFile )

    Parses provided MOL file, extracting molecule information and populating self.

* Molecule.**fGroupSearcher**()

    Scans molecule for functional groups.

#### Mol2D

#### Mol3D


[Atom]: #atom-int-index-array-position-string-element-int-charge--0-
[Bond]: #bond-int-index-atom-bondstart-atom-bondend-int-bondtype--1-int-bonddirection--0-
