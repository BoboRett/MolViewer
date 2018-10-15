# MolViewer
**MolViewer** is a JavaScript library for drawing Chemical Structures. It brings tools to easily convert MOL files and SMILE strings into 2D and 3D representations.

## Resources
* [API Reference](#api)
* [Examples](https://boborett.github.io/MolViewer/)


## API

### Classes
#### Atom( *int* index, *array* position, *string* element[, *int* charge = 0] )

##### Properties
* *int* Atom.**index**

* *obj* Atom.**pos**

    Vector-like: { x, y, z }

* *string* Atom.**element**

* *int* Atom.**charge**

* *array* Atom.**bondedTo**

   Gives an array of connections with Atom, in form of { bond, pairedAtom }

#### Bond( *int* index, *Atom* bondStart, *Atom* bondEnd[, *int* bondType = 1, *int* bondDirection = 0] )

##### Properties
* *int* Bond.**index**

* *[Atom][Atom]* Bond.**start**

* *[Atom][Atom]* Bond.**end**

* *int* Bond.**type**

   * 1: Single
   * 2: Double
   * 3: Triple
   * 4: Aromatic
   * 9: Transition

* *int* Bond.**direction**

   Bond stereochemistry

   * 1: Wedged
   * 6: Hashed

* *bool* Bond.**claimed**

    Used internally during molecular searches, see [Molecule.fGroupSearcher()][Molecule]

#### fGroup( [*[Atom][Atom]* source, *array[[Atom][Atom]]* domain, *array[[Bond][Bond]]* claimed, *string* type ] )

##### Properties

* *string* fGroup.**type**

   Name of functional group, i.e. Primary Amine, Hydroxyl, etc.

* *array[[Atom][Atom]]* fGroup.**domain**

   List of atoms contained within group

* *array[[Bond][Bond]]* fGroup.**claimed**

   List of bonds owned by group

#### Molecule( [*string* molFile] )

##### Properties
* *string* Molecule.**molFile**

   Contains the Molecule's mol file. When setting, the molecule will automatically parse the file and populate itself with atoms, bonds, and functional groups.

* *array[[Atom][Atom]]* Molecule.**atoms**

* *array[[Bond][Bond]]* Molecule.**bonds**

* *array* Molecule.**fGroups**

    Array of molecule's functional groups, in form of {}


##### Methods

* Molecule.**parseMol**( *string* molFile )

    Parses provided MOL file, extracting molecule information and populating self.

* Molecule.**fGroupSearcher**()

    Scans molecule for functional groups.

* Molecule.**get2DFromSMILE**()



#### Mol2D

#### Mol3D


[Atom]: #atom-int-index-array-position-string-element-int-charge--0-
[Bond]: #bond-int-index-atom-bondstart-atom-bondend-int-bondtype--1-int-bonddirection--0-
[Molecule]: #Molecule
