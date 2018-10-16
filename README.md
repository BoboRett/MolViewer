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
        * [Frame Functions](#frame-functions)
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
<div id="mol3D" style="height:500px;width:500px"></div>

<script src="https://d3js.org/d3.v5.js"></script>
<script src="https://cdn.jsdelivr.net/gh/mrDoob/three.js@r97/build/three.min.js"></script>
<script src="https://cdn.jsdelivr.net/gh/BoboRett/MolViewer@v0.52/molViewer.js"></script>
<script src="https://cdn.jsdelivr.net/gh/mrDoob/three.js@r97/examples/js/effects/OutlineEffect.js"></script>
<script src="https://cdn.jsdelivr.net/gh/mrDoob/three.js@r97/examples/js/controls/OrbitControls.js"></script>

<script>
    const Molecule = new MolViewer.Molecule()
    Molecule.get3DFromSMILE( "C1CCC(N)C1" );

    const mol3D = new MolViewer.Mol3D( null, document.getElementById( "mol3D" ) );
    mol3D.init();

    document.addEventListener( "ajaxComplete", () => { mol3D.Molecule = Molecule; mol3D.draw() } );
</script>

```

## Prerequisites

To utilise this entire library, you'll need 5 other JavaScript modules.
* [D3](https://d3js.org/) (Last tested on v5.7.0)

    *Essential* for 2D drawing
* [three.js](https://threejs.org/) (Last tested on r97)

    *Essential* for 3D drawing
    * [three.js - Core](https://github.com/mrdoob/three.js/tree/dev/build)
    * [three.js - OutlineEffect](https://github.com/mrdoob/three.js/tree/dev/examples/js/effects)
    * [three.js - Orbit Controls](https://github.com/mrdoob/three.js/tree/dev/examples/js/controls)

* [OpenChemLib](https://github.com/cheminfo/openchemlib-js) (Last tested on 5.6.0)

    Fairly optional. Just required to convert Smile string into 2D.

* *For page performance stats, you'll need [stats.js](https://github.com/mrdoob/stats.js/). See [Mol3D parameters](#optional-parameters-1)*


## Installing

All you need to do is include the script tag in your HTML.

Github itself can't be used as a CDN (content delivery network), so linking straight to the .js file above won't get you anywhere.

---
The easiest thing is to just yank the script off of a wonderful free CDN that is hooked into Github. Normally, I'd suggest RawGit, but apparently Cryptominers were using it for nefarious deeds, so the next best thing I can find would be [jsDelivr](https://www.jsdelivr.com/).  
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

### ( [ *int* index, *array* position, *string* element, *int* charge = 0 ] )

Container class for Atoms.

### Properties

---

| Property | Type | Description |
| --- | --- | --- |
| index | *int* | Atom index, starts at 0 |
| pos | *obj* | Atom position, in form of `{ x, y, z }` |
| element | *string* | Atom element, i.e. C, N, Cl |
| charge | *int* | Atom parity |
| bondedTo | *array* | Gives an array of connections with Atom, in form of `{ bond: Connecting bond, atom: Paired atom }` |
| object2D | *[SVGSVGElement](https://developer.mozilla.org/en-US/docs/Web/API/SVGSVGElement)* | SVG element that corresponds to Atom |
| object3D | *[Object3D](https://threejs.org/docs/#api/en/core/Object3D)* | Object that corresponds to Atom |
---
---

## Bond

### ( [ *int* index, *Atom* bondStart, *Atom* bondEnd, *int* bondType = 1, *int* bondDirection = 0 ] )

Container class for Bonds.

### Properties

---

| Property | Type | Description |
| --- | --- | --- |
| index | *int* | Bond index, starts at 0 |
| start | *[Atom]* | Starting Atom |
| end | *[Atom]* | Ending Atom |
| type | *int* | Bond type; *1: Single, 2: Double, 3: Triple, 4: Aromatic, 9: Transition* |
| direction | *int* | Bond Stereochemistry; *1: Wedged, 6: Hashed* |
| claimed | *bool* | Used internally during molecular searches, see [Molecule.fGroupSearcher()][Molecule] |
| object2D | *[SVGSVGElement](https://developer.mozilla.org/en-US/docs/Web/API/SVGSVGElement)* | SVG element that corresponds to Bond |
| object3D | *[Object3D](https://threejs.org/docs/#api/en/core/Object3D)* | Object that corresponds to Bond |

---
---

## fGroup

### ( [ *[Atom]* source, *array[[Atom]]* domain, *array[[Bond]]* claimed, *string* type ] )

Container class for functional groups.

### Properties

---

| Property | Type | Description |
| --- | --- | --- |
| type | *string* | Name of functional group, i.e. Primary Amine, Hydroxyl, etc. |
| domain |  *array[[Atom]]* | List of atoms contained within group |
| claimed | *array[[Bond]]* | List of bonds owned by group |
| material | *[Material](https://threejs.org/docs/#api/en/materials/Material)* | Material for functional group encasing mesh. Can be set before drawing or modified afterwards. |

---
---

## Molecule

### ( [ *string* molFile ] )

### Properties

---

| Property | Type | Description |
| --- | --- | --- |
| molFile | *string* | Contains the Molecule's mol file. When setting, the molecule will automatically parse the file and populate itself with atoms, bonds, and functional groups
| atoms | *array[[Atom]]* | List of molecule's atoms
| bonds | *array[[Bond]]* | List of molecule's bonds
| fGroups | *array[[fGroup]]* | List of molecule's functional groups
| bondLength | *int* | Average length of bonds in molecule. Can set this to 'scale' molecule

### Methods

---

| Method | Arguments | Returns | Description |
| --- | --- | --- | --- |
| parseMol | *string* molFile | self | Parses provided MOL file, extracting molecule information and populating self. |
| get2DFromSMILE | *string* smile <br> [*bool* addHydrogrens] | self | Takes **smile** (i.e. "C1CCCC(N)C1" ) and generates and parses its corresponding molfile. Gives only 2D coordinates.<br>Use **addHydrogens** to generate implicit hydrogens with molecule.<br>***Requires OpenChemLib, see [Prerequisites](#prerequisites)*** |
| get3DFromSMILE | *string* smile | void | Takes **smile** (i.e. "C1CCCC(N)C1" ) and generates and parses its corresponding molfile. Gives only 3D coordinates<br>***Synchronously gets conversion via [NIH - National Cancer Institute](https://cactus.nci.nih.gov/translate/). Immediately calling draw will not work, as the fetch will not have completed. You must listen for the "ajaxComplete" event on the document to draw the resulting fetch.*** |
| centre | *none* | void | Centres molecule around 0,0,0 if molfile coordinates are iffy |

---
---

## Mol2D

### ( *[Molecule]* Molecule, *[DOMElement](https://developer.mozilla.org/en-US/docs/Web/API/Element)* Container, *object* dims[, *object* params ] )

### Optional Parameters

| Parameter | Default | Description |
| --- | --- | --- |
| zoomable | true | Enables manipulation of SVG by mouse/touch |
| zoomSmooth | 300 | Duration of zoom transition |
| zoomEase | d3.easeCircleOut | D3 easing function for zoom transition |
| showIndices | false | Show atom indices |
| showHydrogen | true | Show hydrogens |

### Properties

---

| Property | Type | Description |
| --- | --- | --- |
| Molecule | *[Molecule]* | Molecule object to pull information from |
| Container | *[DOMElement](https://developer.mozilla.org/en-US/docs/Web/API/Element)* | Container for SVG |
| dims | *obj* | Dimensions of SVG ViewBox, in form of `{ x, y, width, height }` |
| bondScale | *int* | Another layer of scaling the molecule. See [Molecule] => bondLength |
| stylesheet | *string* | Styles used by 2D molecule. Applied to DOM after .init() |

### Methods

---

| Method | Arguments | Returns | Description |
| --- | --- | --- | --- |
| init | *none* | self | Initialises SVG, adding stylesheet to DOM |
| draw | *none* | self | Renders SVG, adds to .Container |
| genMolecule | *none* | [SVGGElement](https://developer.mozilla.org/en-US/docs/Web/API/SVGGElement) | Generates an SVG Group containing all molecule elements |
| fitToScreen | *none* | void | Fits SVG to current container |

---
---

## Mol3D

### ( *[Molecule]* Molecule, *[DOMElement](https://developer.mozilla.org/en-US/docs/Web/API/Element)* Container[, *object* params ] )

### Optional Parameters

| Parameter | Default | Description |
| --- | --- | --- |
| disableInteractions | false | Disable interactions with 3D canvas |
| showfGroups | true | Show functional groups on molecule |
| showHs | true | Show hydrogens on molecule |
| showStats | false | Show a performance counter in the corner of your page to monitor page stats |
| highlight | true | [Frame function](#frame-functions): Highlights objects upon mouseover |
| autoRotate | false | [Frame function](#frame-functions): Automatically rotates molecule upon axis |
| highlightSync | false | [Frame function](#frame-functions): Highlights corresponding 2D element when hovering a 3D object. **Requires identical atom maps.** |
| mouseoverDispatch | false | [Frame function](#frame-functions): Dispatches events when entering and leaving objects. Carries payload identifying hovered objects. |
| labelTrack | true | [Frame function](#frame-functions): Provides a framework for attaching divs to objects, snapping them during transformations to molecule |

### Properties

---

| Property | Type | Description |
| --- | --- | --- |
| Molecule | *[Molecule]* | Molecule object to pull information from |
| Container | *[DOMElement](https://developer.mozilla.org/en-US/docs/Web/API/Element)* | Container for WebGL canvas |
| Scene | *[Scene](https://threejs.org/docs/#api/en/scenes/Scene)* | Dimensions of SVG ViewBox, in form of `{ x, y, width, height }` |
| molGroup | *[Group](https://threejs.org/docs/#api/en/objects/Group)* | Group of all molecule objects in Scene |
| stylesheet | *string* | Styles used by 3D molecule. Applied to DOM after .init() |
| atomCols | *object(Proxy)* | List of default colours for atoms. See [Colours](#colours) |
| bondCols | *object(Proxy)* | List of default colours for bonds. See [Colours](#colours) |
| groupCols | *object(Proxy)* | List of default colours for functional groups. See [Colours](#colours) |
| labels | *array* | Equivalent to `Mol3D.frameFunctions["labelTrack"].props.labels`<br>See [Frame Functions](#frame-functions) => labelTrack. |


### Methods

---

| Method | Arguments | Returns | Description |
| --- | --- | --- | --- |
| init | *none* | *void* | Initialises scene, adding all required components (lights, cameras---not action, that comes next) |
| draw | *none* | *void* | Draws all objects to scene, adds event listeners for mouse position and window size. Then starts running the scene. |
| play | *none* | *void* | Plays the scene, after being paused. |
| pause | *none* | *void* | Pauses the scene, halting all frameFunctions |
| dispose | *none* | *void* | Removes WebGL canvas from DOM and disposes of scene objects |
| genGroup | *none* | *[Group](https://threejs.org/docs/#api/en/objects/Group)* | Generates a group containing all atoms, bonds, and functional groups
| resetView | *none* | *void* | Resets the camera to fit the molecule into view |
| printHeirarchy | [*[Object3D](https://threejs.org/docs/#api/en/core/Object3D)* object] | *void* | Just a helper function that prints out a tree of all objects and descendants inside of the provided **object**. If none is provided, it uses Mol3D.molGroup |
| setView | *[Object3D](https://threejs.org/docs/#api/en/core/Object3D)* lookAt,<br>*[Vector3](https://threejs.org/docs/#api/en/math/Vector3)* position,<br>*[Vector3](https://threejs.org/docs/#api/en/math/Vector3)* up | *void* | Sets camera view to point at **lookAt** object, and **position**s it, orienting its **up** axis |
| toggleCam | *none* | *void* | Toggles the active camera between perspective and orthographic |


### Frame functions

---
<table>
<thead>
<tr><th>Function Name</th><th>Default Enabled</th><th>Properties</th><th>Default</th><th>Description</th></tr>
</thead>
<tbody>
<tr>
<td>sceneController</td><td>true</td><td>*none*</td><td></td><td>**Required for scene to 'run'. Handle with extreme care.**<br>Renders scene, updates controls and camera frustrums</td>
</tr>
<tr>
<td>highlight</td><td>true</td><td>*none*</td><td></td><td>Highlights atoms as they're hovered</td>
</tr>
<tr>
<td>autoRotate</td><td>false</td><td colspan=2></td><td>Automatically rotates the molecule</td>
</tr>
<tr>
<td rowspan=3 colspan=2></td><td>*string* axes</td><td>x</td><td>Axes upon which to rotate. i.e. "x", "xz", "xyz"</td>
</tr>
<tr>
<td>*<a href="https://threejs.org/docs/#api/en/core/Object3D">Object</a>* target</td><td>Mol3D.molGroup</td><td>Objects that will be rotated about world axis</td>
</tr>
<tr>
<td>*float* speed</td><td>0.01</td><td>Speed at which object rotates (rad/s)</td>
</tr>
<tr>
<td>highlightSync</td><td>false</td><td>*none*</td><td></td><td>Highlights corresponding 2D element when hovering a 3D object.<br>**Requires identical atom maps.**<br>Generating a molecule from get2DFromSMILE and get3DFromSMILE does **not** guarantee identical atom indices, meaning it will highlight the wrong SVG elements. This is really only appropriate when you have imported your own molfiles and can ensure the atoms and bonds are ordered the same way.</td>
</tr>
<tr>
<td>mouseoverDispatch</td><td>false</td><td>*none*</td><td></td><td>Dispatches events to *document* whenever an atom is hovered ("*3DMousein*") or left ("*3DMouseout*"). The event detail contains the **type** of event and an array of **objects** that are under the pointer.</td>
</tr>
<tr>
<td>labelTrack</td><td>false</td><td colspan=2></td><td>Provides a framework for attaching divs to objects, snapping them during transformations to molecule</td>
</tr>
<tr>
<td colspan=2></td><td>*array* labels</td><td>*empty*</td><td>List of label elements with owning 3D **object** attached as data. [~~See Example~~ (coming soon)]()</td>
</tr>
</tbody>
</table>

#### Build your own

Several things rely on calculations each frame. A basic one, that is critical for to work, is simply calling for the next frame to be rendered; the scene and effects need updating. Instead of constantly jacking into the requestAnimationFrame side of things, I just call that once, then iterate each of the functions contained within Mol3D.frameFunctions. The structure is quite specific:

```
Mol3D.frameFunctions = [

    functionName : {

        enabled: boolean - Whether you want it running or not,

        props: obj - Optional, anything internal to the function that needs storing between frames,

        fn: function - The actual function. self is passed as the first argument.

    },

    nextFunctionName : ...,

    ...

]
```

As an example, the autoRotate structure looks as below:

```
...
},
//////Auto rotate//////
autoRotate: {

    enabled: false,

    props: { axes: "y", target: this.molGroup, speed: 0.01 },

    fn: function( self ){

        this.props.target.rotateOnWorldAxis( new THREE.Vector3(
            this.props.axes.includes("x") ? 1 : 0,
            this.props.axes.includes("y") ?  1 : 0,
            this.props.axes.includes("z") ? 1 : 0
        ), this.props.speed )

    }

},

...
```

You can add your own functions to the frameFunction array, just make sure you follow that structure and you'll be grand.


### Colours

---

The colouring of atoms is apparently quite contentious...This library, by default, follows CPK standards, but you can set your own colours with ease.

Every instance of Mol3D contains three collections of colours: **atomCols**, **bondCols**, and **groupCols**. These are all 'Proxies', which mean they contain the usual object structure, but also provide a default option if you try to access something that doesn't exist.

The shortest and easiest to look at is the **bondCols**:

```
this.bondCols = new Proxy( {
    1 : [0.1, 0.1, 0.1],
    2 : [0.5, 0.1, 0.3],
    3 : [0.1, 0.2, 0.5],
    4 : [0.9, 0.5, 0.2],
    9 : [0.3, 0.3, 0.3]
}, {
    get: function( target, name ){

        if( !target.hasOwnProperty( name ) ){ target[name] = [ 0, 0, 0 ] };

        return target[name]

    }
});
```

What this tells you, for instance, is that a triple bond (of Btype = 3), will have a colour in RGB of 0.1, 0.2, 0.5 (![#19327f](https://placehold.it/15/19327f?text=+)). If you tried to get a bond of bType = 11, then the getter will kick in, returning a solid black RGB of 0, 0, 0 (![#000000](https://placehold.it/15/000000?text=+)).

You can modify and set these values as you normally would any value in an object:

```
console.log( bondCols[3] ); // : [0.1, 0.2, 0.5]
bondCols[3] = [1, 1, 1];
console.log( bondCols[3] ); // : [1, 1, 1]

console.log( bondCols[11] ); // : Returns default colour [0, 0, 0]
bondCols[11] = [0.5, 0.5, 0.5 ];
console.log( bondCols[11] ); // : [0.5, 0.5, 0.5]
```

*These colours will only propagate after a .draw(), though.*  
Any dynamic updates to colours can be done by looking at your instance's Colour proxies after drawing:

```
A sample of atomCols after drawing

...
Bk: (3) [0.54, 0.31, 0.89]
Br: (3) [0.65, 0.16, 0.16]
C: (3) [0.56, 0.56, 0.56, material: Ab]
Ca: (3) [0.24, 1, 0]
Cd: (3) [1, 0.85, 0.56]
...
```

When something is drawn for the first time, it will generate a material and attach it to its corresponding value in a Col Proxy (be it Atom, Bond, or Functional Group). That material is referenced in all instances of the object, so modifying the material here *will* propagate dynamically.

If you wanted to completely overwrite the default toon material used, then you could make an entire Material with three.js and attach it to its corresponding Proxy

[Atom]: #atom
[Bond]: #bond
[fGroup]: #fgroup
[Molecule]: #molecule
[Mol2D]: #mol2d
[Mol3D]: #mol3d
