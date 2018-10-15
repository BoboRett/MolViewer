# Intro
**MolViewer** is a JavaScript library for drawing Chemical Structures. It brings tools to easily convert MOL files and SMILE strings into 2D and 3D representations.

# Resources
* [Getting Started](#getting-started)
    * [Quick Start](#quick-start)
    * [Prerequisites](#prerequisites)
    * [Installing](#installing)
* [API Reference](#api)
    * [Atom]
    * [Bond]
    * [fGroup]
    * [Molecule]
    * [Mol2D]
    * [Mol3D]
* [Demo](https://boborett.github.io/MolViewer/)

---
---
---

# Getting Started

These instructions should hopefully get you off the ground with regards to using this library. My style is usually to expose everything (also because I'm lazy), so inspecting the objects might show you a lot of things you're never going to need to touch. However, I designed it as a basic framework for my own tools, so hopefully it's very versatile for anyone else wanting to make something cool with it.

Feel free to send me any examples that could be featured here.

## Quick Start

Without knowing anything, you can get a molecule floating in your webpage with the following:
```
<div id="mol3D"></div>

<script src="https://d3js.org/d3.v5.js"></script>
<script src="https://cdn.jsdelivr.net/gh/mrDoob/three.js@r97/build/three.min.js"></script>
<script src="https://cdn.jsdelivr.net/gh/BoboRett/MolViewer@v0.5/molViewer.js"></script>

<script>
    const Molecule = new MolViewer.Molecule().get3DFromSMILE( "C1CCC(N)C" )

    const Mol3D = new MolViewer.Mol3D( Molecule, document.getElementById( "mol3D" ) ).init().draw()
</script>
```

## Prerequisites

To utilise this entire library, you'll need 3 other JavaScript modules.
* [D3](https://d3js.org/)(Last tested on v5.7.0)

    *Essential* for 2D drawing
* [three.js](https://threejs.org/)(Last tested on r97)

    *Essential* for 3D drawing
* [OpenChemLib](https://github.com/cheminfo/openchemlib-js)(Last tested on 5.6.0)

    Fairly optional. Just required to convert Smile string into 2D.


## Installing

All you need to do is include the script tag in your HTML.

Github itself can't be used as a CDN (content delivery network), so linking straight to the .js file above won't get you anywhere.

---
The easiest thing is to just yank the script off of a wonderful free CDN that is hooked into Github. Normally, I'd suggest RawGit, but apparently CryptoMiners ruined that for everyone, so the next best thing I can find would be [jsDelivr](https://www.jsdelivr.com/).  
A simple template:
```
<script src="https://cdn.jsdelivr.net/gh/BoboRett/MolViewer@INSERTVERSIONNAMEHERE/molViewer.js"></script>
```
---
---
You *can* download the file and host it for your own needs. If you put it in the same folder as your webpage itself, then it's a simple:
```
<script src="molViewer.js"></script>
```
---
Otherwise, you can navigate around relative to your page with a mix of .'s and //'s. For instance, if I had my webpage in a folder, and my scripts in a subfolder, it may look like:

* **MyAwesomeProject**
    * EvenAwesomererWebpage.html
    * **Scripts**
        * molViewer.js

```
<script src="Scripts/molViewer.js"></script>
```
---
Or, if my webpages were kept nested inside a project folder:

* **MyAwesomeProject**
    * **Pages**
        * EvenAwesomererWebpage.html
    * **Scripts**
        * molViewer.js

```
<script src="../Scripts/molViewer.js"></script>
(../ means 'up' a folder to get us to the project's root directory)
```
---

# API

## Atom

### MolViewer.Atom( [*int* index, *array* position, *string* element, *int* charge = 0] )

Container class for Atoms.

### Properties

---

| Property | Type | Description |
| --- | --- | --- |
| index | *int* | Atom index, starts at 0 |
| pos | *obj* | Atom position, in form of { x, y, z } |
| element | *string* | Atom element, i.e. C, N, Cl |
| charge | *int* | Atom parity |
| bondedTo | *array* | Gives an array of connections with Atom, in form of { bond, pairedAtom } |
| object2D | *[SVGSVGElement](https://developer.mozilla.org/en-US/docs/Web/API/SVGSVGElement)* | SVG element that corresponds to Atom |
| object3D | *[three.js Object](https://threejs.org/docs/#api/en/core/Object3D)* | Object that corresponds to Atom |
---
---

## Bond

### MolViewer.Bond( [*int* index, *Atom* bondStart, *Atom* bondEnd, *int* bondType = 1, *int* bondDirection = 0] )

Container class for Bonds.

### Properties

---

### *int* *Bond*.index

### *[Atom]* *Bond*.start

### *[Atom]* *Bond*.end

### *int* *Bond*.type

   * 1: Single
   * 2: Double
   * 3: Triple
   * 4: Aromatic
   * 9: Transition

### *int* *Bond*.direction

   Bond stereochemistry

   * 1: Wedged
   * 6: Hashed

### *bool* *Bond*.claimed

Used internally during molecular searches, see [Molecule.fGroupSearcher()][Molecule]

### *[SVGSVGElement](https://developer.mozilla.org/en-US/docs/Web/API/SVGSVGElement)* *Bond*.object2D

SVG element that corresponds to Bond

### *[three.js Object](https://threejs.org/docs/#api/en/core/Object3D)* *Bond*.object3D

Object that corresponds to Bond

---
---

## fGroup

### MolViewer.fGroup( [*[Atom]* source, *array[[Atom]]* domain, *array[[Bond]]* claimed, *string* type ] )

Container class for functional groups.

### Properties
---

### *string* *fGroup*.type

   Name of functional group, i.e. Primary Amine, Hydroxyl, etc.

### *array[[Atom]]* *fGroup*.domain

   List of atoms contained within group

### *array[[Bond]]* *fGroup*.claimed

   List of bonds owned by group

---
---

## Molecule

### MolViewer.Molecule( [*string* molFile] )

### Properties
---

### *string* *Molecule*.molFile

Contains the Molecule's mol file. When setting, the molecule will automatically parse the file and populate itself with atoms, bonds, and functional groups.

### *array[[Atom]]* *Molecule*.atoms

### *array[[Bond]]* *Molecule*.bonds

### *array[[fGroup]]* *Molecule*.fGroups

Array of molecule's functional groups, in form of {}

### *int* *Molecule*.bondLength

Average length of bonds in molecule. Can change this to 'scale' molecule.


### Methods

---

### *Molecule*.parseMol( *string* molFile )

Parses provided MOL file, extracting molecule information and populating self.  
Returns array of found [[Atom]s, [Bond]s]

### *Molecule*.fGroupSearcher()

Scans molecule for functional groups.  
Returns array of found [fGroup]s

### *Molecule*.get2DFromSMILE( *string* smile[, *bool* addHydrogens ] )

Takes **smile** (i.e. "C1CCCC(N)C1" ) and generates and parses its corresponding molfile. Gives only 2D coordinates.
Use **addHydrogens** to generate implicit hydrogens with molecule.  
***Requires OpenChemLib, see [Prerequisites](#prerequisites)***  

### *Molecule*.get3DFromSMILE( *string* smile )

Takes **smile** (i.e. "C1CCCC(N)C1" ) and generates and parses its corresponding molfile. Gives only 3D coordinates  
    ***Gets conversion via [NIH - National Cancer Institute](https://cactus.nci.nih.gov/translate/)***  

### *Molecule*.centre()

  Centres molecule around 0,0,0 if molfile coordinates are iffy

---
---

## Mol2D

### MolViewer.Mol2D( *Molecule* molecule, *[DOMElement](https://developer.mozilla.org/en-US/docs/Web/API/Element)* container, *object-like* dims[, *object* params ] )

---
---

## Mol3D

### MolViewer.Mol3D( *Molecule* molecule, *[DOMElement](https://developer.mozilla.org/en-US/docs/Web/API/Element)* container[, *object* params ] )


[Atom]: #atom
[Bond]: #bond
[fGroup]: #fgroup
[Molecule]: #molecule
[Mol2D]: #mol2d
[Mol3D]: #mol3d
