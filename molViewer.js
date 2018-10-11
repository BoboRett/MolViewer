//Molecule Viewer. Author: Jamie Ridley
(function( global, factory ){

	factory((global.MolViewer = {}))

}(this, function( exports ){

	/////Atom object/////
	function Atom( index, position, element, charge ){

		this.index = index;
		this.pos = position;
		this.element = element;
		this.charge = charge;
		this.bondedTo = [];

	}

	/////Bond object/////
	function Bond( index, bondStart, bondEnd, bondType, bondDirection ){

		this.index = index;
		this.start = bondStart;
		this.end = bondEnd;
		this.type = bondType;
		this.direction = bondDirection;
		this.claimed = false;

	}

	/////Molecule object/////
	function Molecule( molData ){

		if( molData ){

			this.molFile = molData;

		} else{

			this._molFile = null;
			this.atoms = [];
			this.bonds = [];
			this.fGroups = [];

		}

	};

	/////Molecule methods/////
	Object.assign( Molecule.prototype, {

		parseMol: function( molData ){

			let mol = {};
			molData = molData.split( "\n" ).map( el => el.trim() ).join( "\n")
			const [numatoms,numbonds] = molData.split( "\n" )[3].match( /.{1,3}/g ).slice( 0, 2 ).map( el => parseInt( el ) );

			mol.atoms = molData.split( "\n" )
							.slice( 4, 4 + numatoms )
							.map( function( el, i ){

								el = " ".repeat( 69 - el.length ) + el;
								const line = el.slice( 0, 30 ).match( /.{1,10}/g ).concat( el.slice( 30 ).match( /.{1,3}/g ) );
								return new Atom(
										i,
										[parseFloat( line[0] ),parseFloat( line[1] ), parseFloat( line[2] )],
										line[3].trim(),
										0
									)

							});
			mol.bonds = molData.split( "\n" )
							.slice( 4 + numatoms, 4 + numatoms + numbonds )
							.map( function( el, i ){

								el = " ".repeat( 21 - el.length ) + el;
								const line = el.match( /.{1,3}/g );
								return new Bond(
										i,
										mol.atoms[line[0] - 1],
										mol.atoms[line[1] - 1],
										parseInt( line[2] ),
										parseInt( line[3] ),
									)

							});

			//////CHARGES//////
			molData.split( "\n" ).forEach( function( el, i ){
				const line = el.match( /\S+/g )
				if( line !== null ){

					if( line[1] === "CHG" ){

						for( var i = 3; i < line.length; i = i + 2 ){

							mol.atoms[ +line[i] - 1 ].charge = +line[i + 1]

						}

					}

				}
			})

			//////ATTACH BONDS TO ATOMS//////
			mol.atoms.forEach( function( atom, i ){
				mol.bonds
					.filter( bond => bond.start.index === i || bond.end.index === i )
					.forEach( function( bond ){
						var from = i;
						var to = ( bond.start.index === i ? bond.end.index : bond.start.index );
						atom.bondedTo.push( {el: mol.atoms[to], bond: bond} );
					})
			})

			return [mol.atoms, mol.bonds];
		},

		fGroupSearcher: function( mol ){

			let fGroups = []

			let subStructures = [
								 {type: "Carboxylic Acid", root: "C", bonds:[{el: "O", btype: "2"}, {el: "O", btype: "1", bondedTo:[{el: "H", btype: "1"}]}]},
								 {type: "Ester", root: "C", bonds:[{el: "O", btype: "2"}, {el: "O", btype: "1", bondedTo:[{el: "R", btype: "1"}]}]},
								 {type: "Amide", root: "C", bonds:[{el: "O", btype: "2"}, {el: "N", btype: "1", bondedTo:[{el: "R", btype: "1"}, {el: "R", btype: "1"}]}]},

								 {type: "Acyl Halide", root: "C", bonds:[{el: "O", btype: "2"}, {el: "X", btype: "1"}]},

								 {type: "Aldehyde", root: "C", bonds:[{el: "O", btype: "2"}, {el: "H", btype: "1"}]},
								 {type: "Ketone", root: "C", bonds:[{el: "O", btype: "2"}]},

								 {type: "Primary Amine", root: "N", bonds:[{el: "H", btype: "1"},{el: "H", btype: "1"}]},
								 {type: "Secondary Amine", root: "N", bonds:[{el: "H", btype: "1"},{el: "R", btype: "1"}]},
								 {type: "Tertiary Amine", root: "N", bonds:[{el: "R", btype: "1"},{el: "R", btype: "1"}]},

								 {type: "Nitro", root: "N", bonds:[{el: "O", btype: "2"},{el: "O", btype: "1"}]},
								 {type: "Alcohol", root: "O", bonds:[{el: "H", btype: "1"}]},
								 {type: "Halo", root: "R", bonds:[{el: "X", btype: "1"}]},
								 {type: "Nitrile", root: "C", bonds:[{el: "N", btype: "3"}]},

								 {type: "Acetal", root: "C", bonds:[{el: "O", btype: "1"}, {el: "O", btype: "1"}]},
								 {type: "Ether", root: "O", bonds:[{el: "R", btype: "1"},{el: "R", btype: "1"}]},
								];

			mol.atoms.filter( atom => atom.element != "H" ).forEach( function( atom ){
				let scanning = true

				while( scanning ){
					scanning = false

					for( ss of subStructures.filter( ss => ( ss.root === "R" ? true : ss.root === atom.element ) ) ) { //Filter out non-matching root atoms

						let foundStructures = inBonds( atom, ss.bonds, atom.index, [] );

						if( foundStructures.hasOwnProperty( "source" ) ){
							foundStructures.type = ss.type;
							fGroups = fGroups.concat( foundStructures );
							foundStructures.claimed.map( el => {mol.bonds[el.index].claimed = true} )
							scanning = true;
							break;
						}

					}
				}

			})

			//////SEARCH BONDS OF ATOM//////
			function inBonds( source, subStruct, rootIndex, domain ){

				let claimed = [];

				subStruct.forEach( function( ss ){

					ss.found = false;
					source.bondedTo.filter( bond => bond.el.index !== rootIndex && !domain.map( dom => dom.index ).includes( bond.el.index ) && bond.bond.claimed === false ).forEach( function( bond ){
						if( !ss.found && bond.bond.type === ss.btype ){
							if( ( ss.el === "R" ? true : ( ss.el === "X" ? ["Cl", "Br", "I", "F"].includes( bond.el.element ) : bond.el.element === ss.el ) ) ){

								domain.push( bond.el );
								claimed.push( bond.bond );

								if( ss.hasOwnProperty( "bondedTo" ) ){

									//////Recursive search//////
									const deepSearch = inBonds( bond.el , ss.bondedTo, rootIndex, domain );

									if( deepSearch.hasOwnProperty( "source" ) ){
										claimed = claimed.concat( deepSearch.claimed );
										ss.found = true;
									}

								} else{
									ss.found = true;
								}

							}
						}
					})

				})

				if( subStruct.filter( el => !el.found ).length < 1 ){
					results = {source: source, domain: domain, claimed: claimed};
				} else{
					results = []
				}

				return results

			}

			return fGroups

		},

		get2DFromSMILE: function( smile, addH ){

			let molecule = OCL.Molecule.fromSmiles( smile );
			addH && molecule.addImplicitHydrogens();
			this.molfile = this.parse( molecule.toMolfile() );

			return this.Mol2D

		},

		get3DFromSMILE: function( smile ){

			this.ajaxRunning = true;
			mol = this;
			const event = new Event( "ajaxComplete" );

			d3.request( "https://cactus.nci.nih.gov/chemical/structure/" + smile.replace( /\#/g, "%23" ).replace( /\[/g, "%5B" ).replace( /\]/g, "%5D" ) + "/file/xml?format=sdf&get3d=true")
				.get( function( err, d ){
					mol.molfile = d.responseXML.children[0].children[0].children[0].innerHTML;
					mol.ajaxRunning = false
					document.dispatchEvent( event )
				})

		}



	})

	/////Molecule properties/////
	Object.defineProperties( Molecule.prototype , {

		"molFile": {

			get: function () {

				return this._molFile;

			},

			set: function ( value ) {

				[this.atoms, this.bonds] = this.parseMol( value );
				this.fGroups = this.fGroupSearcher( this );
				this._molFile = value;

			}

		}

	})

	/////2D canvas objects/////
	function Mol2D( molecule, container, dims, params ){

		params = params || {};

		const self = this;

		this.DOM = container;
		this.dims = dims;
		this.molecule = molecule;

		this.zoomable = params.zoomable || true;
		this.showIndices = params.showIndices || false;

		this.zoomFunc = d3.zoom().on( "zoom", function(){

			self.root.transition()
				.duration( 300 )
				.ease( d3.easeCircleOut )
				.attr( "transform", d3.event.transform )

		})

		this.stylesheet = `
			#view2d {
				text-anchor: middle;
				font-family: sans-serif;
				font-size: 16px;
				cursor: all-scroll;
			}

			.bond > line, .bond_dbl > line, .bond_hash > line, .bond_trp > line {
				stroke: black;
				stroke-width: 1px;
				pointer-events: none;
				stroke-linecap: round;
			}

			.highlight, .highlight_hover {
				stroke: black;
				fill: yellow;
			}

			.highlight {
				stroke-width: 0;
				fill-opacity: 0;
			}

			.highlight_hover {
				stroke-width: 0.5px;
				fill-opacity: 0.5;
			}

			.atomFocus {
				stroke-width: 0.5px;
				stroke: black;
				fill: orange;
			}

			.wedge {
				fill: black;
			}

			.label_O{
				fill: red;
			}
			.label_N{
				fill: blue;
			}

			.atomind > text {
				font-size: 8px;
			}

			.atomind > rect{
				fill: white;
				fill-opacity: 0.7;
			}

			.charge > text{
				text-anchor: middle;
				pointer-events: none;
				font-family: sans-serif;
				font-size: 8px;
				alignment-baseline: middle;
			}

			.charge > circle{
				fill: white;
				fill-opacity: 0.7;
				stroke: black;
				stroke-width: 0.5px;
			}
		`,

		d3.select( "#molViewer2DCSS" ).empty() && d3.select( "head" ).append( "style" ).attr( "id", "molViewer2DCSS" ).html( this.stylesheet );

	};

	/////2D canvas methods/////
	Object.assign( Mol2D.prototype, {

		draw: function(){

			const atoms = this.molecule.atoms;
			const bonds = this.molecule.bonds;

			const svg = this.DOM.append( "svg" ).attr( "viewBox", this.dims.join( " " ) ).attr( "id", "view2d" );
			const root = svg.append( "g" )
					.attr( "id", "rootframe" )
					.attr( "transform", null )
					.html( "" );

			this.svg = svg;
			this.root = root;

			//////ATOMS//////
			const atomsroot = root.append( "g" ).attr( "class", "atoms" )
			atoms.forEach( function( atom, i ){

				const atomGrp = atomsroot.append( "g" )
									.attr( "class", "atom_" + atom.element )
									.attr( "id", i )
									.datum( atom )

				atom.element != "H" && atomGrp.append( "g" ).attr( "class", "hydrogens" );

				const highlightCircle = atomGrp.append( "circle" )
									.attr( "class", "highlight" )
									.attr( "id", "highlight_" + i )
									.attr( "cx", atom.pos[1] )
									.attr( "cy", atom.pos[2] )
									.attr( "r", 12 )

				const txt = atomGrp.append( "text" )
									.attr( "class", "label_" + atom.element )
									.attr( "id", "label_" + i )
									.attr( "x", atom.pos[1] )
									.attr( "y", atom.pos[2] + 6 )
									.text( atom.element !== "C" || atom.charge ? atom.element : "" )

				const ind = atomGrp.append( "g" )
									.attr( "class", "atomind" )
									.attr( "id", "atomind_" + i )
									.attr( "display", this.showIndices ? null : "none" )

				const indBBox = ind.append( "text" )
									.attr( "x", atom.pos[1] - txt.node().getBBox().width/2 - i.toString().split("").length * 2 )
									.attr( "y", atom.pos[2] - 5 )
									.text( i )
									.node().getBBox()

				ind.append( "rect" )
					.attr( "x", indBBox.x )
					.attr( "width", indBBox.width )
					.attr( "y", indBBox.y )
					.attr( "height", indBBox.height )
					.lower();

				if( atom.charge ){

					const chg = atomGrp.append( "g" )
									.attr( "class", "charge" );
					const chgTxt = chg.append( "text" )
										.attr( "class", "atomcharge" )
										.attr( "id", "atomchg_" + i )
										.attr( "x", atom.pos[1] + txt.node().getBBox().width/2 + atom.charge.toString().length * 2 - 2 )
										.attr( "y", atom.pos[2] - 8 )
										.text( atom.charge === -1 ? "-" : ( atom.charge === 1 ? "+" : atom.charge ) )

					chg.append( "circle" )
						.attr( "cx", chgTxt.attr( "x" ) )
						.attr( "cy", chgTxt.attr( "y") )
						.attr( "r", chgTxt.node().getBBox().width/2 + 1 )
						.lower();

				}

				atom.HTML = atomGrp.node();
				atom.highlightCircle = highlightCircle.node()

			})

			//////BONDS//////
			const bondsroot = root.append( "g" ).attr( "class", "bonds" ).lower();
			const labelOffset = 8;

			bonds.forEach( function( bond, j ){

				const tmp = bondsroot.append( "g" ).attr( "class", "bond_" + bond.type ).datum( bond );
				const theta = Math.atan2( bond.end.pos[2] - bond.start.pos[2], bond.end.pos[1] - bond.start.pos[1] );
				const length = Math.hypot( bond.end.pos[1] - bond.start.pos[1], bond.end.pos[2] - bond.start.pos[2] );

				const highlight = tmp.append( "rect" )
					.attr( "x", -8 ).attr( "y", -7.5 )
					.attr( "rx", 7.5 ).attr( "ry", 7.5 )
					.attr( "width", length + 2*8 ).attr( "height", 15)
					.attr("transform", "translate(" + bond.start.pos[1] + "," + bond.start.pos[2] + ")rotate(" + theta*180/Math.PI + ")" )
					.attr( "class", "highlight" )
					.attr( "id", "highlight_" + bond.start.index + "_" + bond.end.index );

				drawBond( tmp, bond, theta, length );

				if( bond.start.element === "H" || bond.end.element === "H" ){

					svg.node().getElementById( bond.start.element === "H" ? bond.end.index : bond.start.index ).getElementsByClassName("hydrogens")[0].appendChild( tmp.node() )
					svg.node().getElementById( bond.start.element === "H" ? bond.end.index : bond.start.index ).getElementsByClassName("hydrogens")[0].appendChild( svg.node().getElementById( ( bond.start.element === "H" ? bond.start.index : bond.end.index ) ) )

				}

				bond.HTML = tmp.node();
				bond.highlight = highlight.node();

			});

			this.zoomable && this.fitToScreen();

			function drawBond( parent, bond, angle, length ){

				let tmp;

				const coords = [bond.start.pos[1] + ( bond.start.element !== "C" || bond.start.charge ? bond.start.element.length * 8 * Math.cos( angle ) : 0 ),
								bond.end.pos[1] - ( bond.end.element !== "C" || bond.end.charge ? bond.end.element.length * 8 * Math.cos( angle ) : 0 ),
								bond.start.pos[2] + ( bond.start.element !== "C" || bond.start.charge ? bond.start.element.length * 8 * Math.sin( angle ) : 0 ),
								bond.end.pos[2] - ( bond.end.element !== "C" || bond.end.charge ? bond.end.element.length * 8 * Math.sin( angle ) : 0 )]

				const placeholderLine = tmp.append( "line" )
						.attr( "class", "bondline")
						.attr( "x1", coords[0] )
						.attr( "x2", coords[1] )
						.attr( "y1", coords[2] )
						.attr( "y2", coords[3] )

				switch( bond.type ){

					case "1": //single bond

						switch( bond.direction ){

							case 0: //normal bond

								parent.attr( "class", "bond" );
								break;

							case 1: //wedge bond

								tmp = parent.attr( "class", "bond_wedge" ).append( "polygon" )
									.attr( "points", coords[0] + "," + coords[2] + " " +
										 ( coords[1] + 3*Math.cos( angle + Math.PI/2 ) ).toFixed( 2 ) + "," + ( coords[3] + 3*Math.sin( angle + Math.PI/2 ) ).toFixed( 2 ) + " " +
										 ( coords[1] - 3*Math.cos( angle + Math.PI/2 ) ).toFixed( 2 ) + "," + ( coords[3] - 3*Math.sin( angle + Math.PI/2 ) ).toFixed( 2 ) );
								placeholderLine.remove()
								break;

							case 6: //hash bond

								tmp = parent.attr( "class", "bond_hash" );
								const point = [0, 0];
								for( let i = 0; Math.hypot( ...point ) < length ; i++ ){
									tmp.append( "line" ).attr( "class", "bond" )
										.attr( "x1", ( coords[0] + point[0] + Math.hypot( ...point )/length * 3*Math.cos( angle + Math.PI/2 ) ).toFixed( 2 ) )
										.attr( "x2", ( coords[0] + point[0] - Math.hypot( ...point )/length * 3*Math.cos( angle + Math.PI/2 ) ).toFixed( 2 ) )
										.attr( "y1", ( coords[2] + point[1] + Math.hypot( ...point )/length * 3*Math.sin( angle + Math.PI/2 ) ).toFixed( 2 ) )
										.attr( "y2", ( coords[2] + point[1] - Math.hypot( ...point )/length * 3*Math.sin( angle + Math.PI/2 ) ).toFixed( 2 ) );
									point = [( i + 1 ) * 3 * Math.cos( angle ), ( i + 1 ) * 3 * Math.sin( angle )];
								};
								placeholderLine.remove()
								break;

						}

						break;

					case "2": //double bond

						parent.attr( "class", "bond_dbl" );
						[-1, 1].forEach( el => {

							parent.node().appendChild( d3.select( placeholderLine.node().cloneNode() )
								.attr( "x1", coords[0] + el*2.5*Math.cos( angle + Math.PI/2 ) )
								.attr( "x2", coords[1] + el*2.5*Math.cos( angle + Math.PI/2 ) )
								.attr( "y1", coords[2] + el*2.5*Math.sin( angle + Math.PI/2 ) )
								.attr( "y2", coords[3] + el*2.5*Math.sin( angle + Math.PI/2 ) )
								.node() );

						})
						placeholderLine.remove()
						break;

					case "3": //triple bond

						parent.attr( "class", "bond_trp" );
						[-1, 1].forEach( el => {

							parent.node().appendChild( d3.select( placeholderLine.node().cloneNode() )
								.attr( "x1", coords[0] + el*2.5*Math.cos( angle + Math.PI/2 ) )
								.attr( "x2", coords[1] + el*2.5*Math.cos( angle + Math.PI/2 ) )
								.attr( "y1", coords[2] + el*2.5*Math.sin( angle + Math.PI/2 ) )
								.attr( "y2", coords[3] + el*2.5*Math.sin( angle + Math.PI/2 ) )
								.node() );

						})
						break;

					case "9": //aromatic bond

						placeholderLine.attr("class","bond")
						break;

				}

			}

			return root
		},

		fitToScreen: function(){
			this.root.attr( "transform", null ) ;
			const viewBox = this.root.node().parentNode.getBoundingClientRect();
			const rootBox = this.root.node().getBoundingClientRect();

			const zoom = viewBox.width/rootBox.width < viewBox.height/rootBox.height ? viewBox.width/rootBox.width : viewBox.height/rootBox.height

			this.svg.call( this.zoomFunc.transform, d3.zoomIdentity.translate( viewBox.width/2 + ( - ( rootBox.left - viewBox.left ) - rootBox.width/2 )*zoom, viewBox.height/2 + ( - ( rootBox.top - viewBox.top ) - rootBox.height/2 )*zoom ).scale( zoom ) )
			this.svg.call( this.zoomFunc );
		},

		showH: function( showH ){

			self.Mol2D.root.selectAll( ".hydrogens, .hydrogens > *" ).each( function(){

				d3.select( this ).attr( "display", showH ? "all" : "none" );

			});

		},

		showLabels: function( show ){

			d3.selectAll( ".atomind" ).attr( "display", show ? null : "none" );

		},

		onWindowResize: function(){
			if( this.svg ){
				this.svg.node().viewBox.baseVal.height = this.svg.node().viewBox.baseVal.width / ( this.DOM.style( "width" ).slice( 0, -2 ) / this.DOM.style( "height" ).slice( 0, -2 ) );
			} else{
				console.warn( "No SVG Element to resize. Call .draw() first!" );
			}

		},

	})

	/////2D canvas properties/////
	Object.defineProperties( Mol2D.prototype, {

		"zoomable": {

			get: function () {

				return this._zoomable;

			},

			set: function ( value ) {

				this._zoomable = value;

			}

		},

		"showIndices": {

			get: function () {

				return this._showIndices;

			},

			set: function ( value ) {

				this._showIndices = value;

			}

		}

	})

	/////3D canvas object/////
	function Mol3D( molecule, container, params ){

		params = params || {};

		const self = this;

		this._initialised = false;
		this._FOV = 70;

		this.Container    = container || null;
		this.Molecule     = molecule;
		this.scene        = new THREE.Scene();
		this.mouse        = new THREE.Vector2();
		this.molGroup     = new THREE.Group();
		this.animID;

		this.stylesheet = `
			.view3D {
				text-anchor: middle;
				font-family: sans-serif;
				font-size: 16px;
			}
		`

		d3.select( "#molViewer3DCSS" ).empty() && d3.select( "head" ).append( "style" ).attr( "id", "molViewer3DCSS" ).html( this.stylesheet );

		this.frameFunctions = {
			//////Scene controller (ONLY TOUCH IF YOU KNOW WHAT YOU'RE DOING)//////
			sceneController: {enabled: true, props: {}, fn: function( self ){

					self.animID = requestAnimationFrame( self._animate );
					self.controls.update();

					self._camOrtho.position.copy( new THREE.Vector3().copy( self._camPersp.position ) );
					self._camOrtho.quaternion.copy( self._camPersp.quaternion );
					let frustum = self._frustum;

					self._camOrtho.left   = -frustum[0] / 2;
					self._camOrtho.right  =  frustum[0] / 2;
					self._camOrtho.top    =  frustum[1] / 2;
					self._camOrtho.bottom = -frustum[1] / 2;

					self._camOrtho.updateProjectionMatrix();

					self.showStats && self.stats.update()

				}
			},
			//////Highlight elements on hover//////
			highlight: {enabled: false, props: { _intersected: null }, fn: function( self ){

					var raycaster = new THREE.Raycaster();
					raycaster.setFromCamera( self.mouse, self.camActive );

					const intersects = raycaster.intersectObjects( self.molGroup.children, true );

					if( self.mouse.x > -1 && self.mouse.x < 1 && self.mouse.y > -1 && self.mouse.y < 1 && intersects.length > 0){

						if( this.props._intersected ){

							this.props._intersected.material = this.props._intersected.currentMat;
							this.props._intersected.material.dispose();

						};

						if( intersects[0].object instanceof THREE.Mesh ){

							this.props._intersected = intersects[0].object;
							this.props._intersected.currentMat = this.props._intersected.material;
							const mat = this.props._intersected.currentMat.clone()
							mat.emissive.setHex( 0xff0000 )
							this.props._intersected.material = mat;

						}
						else{

							this.props._intersected = null;

						}

					} else {

						if( this.props._intersected ){

							this.props._intersected.material = this.props._intersected.currentMat;
							this.props._intersected.material.dispose();

						};

					}

				}
			},
			//////Auto rotate//////
			autoRotate: {enabled: false, props: { axes: "y", target: this.molGroup, speed: 0.01 }, fn: function( self ){
					this.props.target.rotateOnWorldAxis( new THREE.Vector3(
						this.props.axes.includes("x") ? 1 : 0,
						this.props.axes.includes("y") ?  1 : 0,
						this.props.axes.includes("z") ? 1 : 0
					), this.props.speed )
				}
			},
			//////Highlight svg elements from 3D//////
			highlightSync: {enabled: false, props: {}, fn: function( self ){
					var raycaster = new THREE.Raycaster();
					raycaster.setFromCamera( self.mouse, self.camActive );

					var intersects = raycaster.intersectObjects( self.molGroup.children, true );

					if( self.mouse.x > -1 && self.mouse.x < 1 && self.mouse.y > -1 && self.mouse.y < 1 && intersects.length > 0){

						self.hovered && self.hovered.attr( "class", "highlight" );

						var highlighted = intersects[0].object;
						switch( highlighted.userData.type ){

							case "fGroup":

								self.hovered = d3.selectAll( "#highlight_" + highlighted.userData.source.source.index + ", " + highlighted.userData.source.domain.map( el =>  "#highlight_" + el.index ).join( ", " ) )
								break;

							case "atom":
							case "bond":
								self.hovered = highlighted.name.toString().includes("_") ?
									d3.select( "#highlight_" + highlighted.name + ", #highlight_" + highlighted.name.toString().split("").reverse().join("") ) :
									d3.select( "#highlight_" + highlighted.name )
								break;

						}

						self.hovered && self.hovered.attr( "class", "highlight_hover" );

					} else {

						self.hovered && self.hovered.attr( "class", "highlight" )

					}
				}
			},
			//////Dispatch object mouseover events to "3DMouseover"//////
			mouseoverDispatch: {enabled: false, props: { XRay: false, _hovered: null }, fn: function( self ){

					if( !this.props.XRay ){

						if( !self._DOM.matches( ":hover" ) ){ return }

					}

					var raycaster = new THREE.Raycaster();
					raycaster.setFromCamera( self.mouse, self.camActive );

					var intersects = raycaster.intersectObjects( self.molGroup.children, true );

					if( this.props._hovered !== intersects.length && this.props._hovered > 0 ){

						document.dispatchEvent( new CustomEvent( "3DMouseout", {"detail": {"type": "mouseout", "objects":intersects } } ) );

					}

					if( intersects.length > 0 ){
						if( this.props._hovered < intersects.length ){
							document.dispatchEvent( new CustomEvent( "3DMousein", {"detail": {"type": "mouseover", "objects":intersects } } ) );
							this.props._hovered = intersects.length;
						} else if( this.props._hovered > intersects.length ){
							document.dispatchEvent( new CustomEvent( "3DMousein", {"detail": {"type": "mouseout", "objects":intersects } } ) );
							this.props._hovered = intersects.length;
						}
					} else{
						this.props._hovered = 0;
					}
				}
			},
			//////Attach labels to atoms in 3d//////
			labelTrack: {enabled: false, props: { labels: [] }, fn: function( self ){
					var widthHalf = self._DOM.getBoundingClientRect().width/2;
					var heightHalf = self._DOM.getBoundingClientRect().height/2;
					for( label of this.props.labels ){
						var vector = new THREE.Vector3().setFromMatrixPosition( label.datum().object.matrixWorld );
						vector.project( self.camActive );

						label.style("right", ( - ( vector.x * widthHalf ) + widthHalf ) + "px" )
						label.style("bottom", ( ( vector.y * heightHalf ) + heightHalf ) + "px" )
					}
				}
			}
		}

		this.highlight           = params.highlight || true;
		this.autoRotate          = params.autoRotate || false;
		this.highlightSync       = params.highlightSync || false;
		this.mouseoverDispatch   = params.mouseoverDispatch || false;
		this.labelTrack          = params.labelTrack || true;

		this.disableInteractions = params.disableInteractions || false;
		this.showfGroups         = params.showfGroups || true;
		this.showHs              = params.showHs || true;
		this.showStats           = params.showStats || false;

		//////CPK Colours//////
		{
			this.atomCols = new Proxy( {
				H: 		[1,1,1],
				He: 	[0.85,1,1],
				Li: 	[0.8,0.5,1],
				Be: 	[0.76,1,0],
				B: 		[1,0.71,0.71],
				C: 		[0.56,0.56,0.56],
				N: 		[0.19,0.31,0.97],
				O: 		[1,0.05,0.05],
				F: 		[0.56,0.88,0.31],
				Ne: 	[0.7,0.89,0.96],
				Na: 	[0.67,0.36,0.95],
				Mg: 	[0.54,1,0],
				Al: 	[0.75,0.65,0.65],
				Si: 	[0.94,0.78,0.63],
				P: 		[1,0.5,0],
				S: 		[1,1,0.19],
				Cl: 	[0.12,0.94,0.12],
				Ar: 	[0.5,0.82,0.89],
				K: 		[0.56,0.25,0.83],
				Ca: 	[0.24,1,0],
				Sc: 	[0.9,0.9,0.9],
				Ti: 	[0.75,0.76,0.78],
				V: 		[0.65,0.65,0.67],
				Cr: 	[0.54,0.6,0.78],
				Mn: 	[0.61,0.48,0.78],
				Fe: 	[0.88,0.4,0.2],
				Co: 	[0.94,0.56,0.63],
				Ni: 	[0.31,0.82,0.31],
				Cu: 	[0.78,0.5,0.2],
				Zn: 	[0.49,0.5,0.69],
				Ga: 	[0.76,0.56,0.56],
				Ge: 	[0.4,0.56,0.56],
				As: 	[0.74,0.5,0.89],
				Se: 	[1,0.63,0],
				Br: 	[0.65,0.16,0.16],
				Kr: 	[0.36,0.72,0.82],
				Rb: 	[0.44,0.18,0.69],
				Sr: 	[0,1,0],
				Y: 		[0.58,1,1],
				Zr: 	[0.58,0.88,0.88],
				Nb: 	[0.45,0.76,0.79],
				Mo: 	[0.33,0.71,0.71],
				Tc: 	[0.23,0.62,0.62],
				Ru: 	[0.14,0.56,0.56],
				Rh: 	[0.04,0.49,0.55],
				Pd: 	[0,0.41,0.52],
				Ag: 	[0.75,0.75,0.75],
				Cd: 	[1,0.85,0.56],
				In: 	[0.65,0.46,0.45],
				Sn: 	[0.4,0.5,0.5],
				Sb: 	[0.62,0.39,0.71],
				Te: 	[0.83,0.48,0],
				I: 		[0.58,0,0.58],
				Xe: 	[0.26,0.62,0.69],
				Cs: 	[0.34,0.09,0.56],
				Ba: 	[0,0.79,0],
				La: 	[0.44,0.83,1],
				Ce: 	[1,1,0.78],
				Pr: 	[0.85,1,0.78],
				Nd: 	[0.78,1,0.78],
				Pm: 	[0.64,1,0.78],
				Sm: 	[0.56,1,0.78],
				Eu: 	[0.38,1,0.78],
				Gd: 	[0.27,1,0.78],
				Tb: 	[0.19,1,0.78],
				Dy: 	[0.12,1,0.78],
				Ho: 	[0,1,0.61],
				Er: 	[0,0.9,0.46],
				Tm: 	[0,0.83,0.32],
				Yb: 	[0,0.75,0.22],
				Lu: 	[0,0.67,0.14],
				Hf: 	[0.3,0.76,1],
				Ta: 	[0.3,0.65,1],
				W: 		[0.13,0.58,0.84],
				Re: 	[0.15,0.49,0.67],
				Os: 	[0.15,0.4,0.59],
				Ir: 	[0.09,0.33,0.53],
				Pt: 	[0.82,0.82,0.88],
				Au: 	[1,0.82,0.14],
				Hg: 	[0.72,0.72,0.82],
				Tl: 	[0.65,0.33,0.3],
				Pb: 	[0.34,0.35,0.38],
				Bi: 	[0.62,0.31,0.71],
				Po: 	[0.67,0.36,0],
				At: 	[0.46,0.31,0.27],
				Rn: 	[0.26,0.51,0.59],
				Fr: 	[0.26,0,0.4],
				Ra: 	[0,0.49,0],
				Ac: 	[0.44,0.67,0.98],
				Th: 	[0,0.73,1],
				Pa: 	[0,0.63,1],
				U: 		[0,0.56,1],
				Np: 	[0,0.5,1],
				Pu: 	[0,0.42,1],
				Am: 	[0.33,0.36,0.95],
				Cm: 	[0.47,0.36,0.89],
				Bk: 	[0.54,0.31,0.89],
				Cf: 	[0.63,0.21,0.83],
				Es: 	[0.7,0.12,0.83],
				Fm: 	[0.7,0.12,0.73],
				Md: 	[0.7,0.05,0.65],
				No: 	[0.74,0.05,0.53],
				Lr: 	[0.78,0,0.4],
				Rf: 	[0.8,0,0.35],
				Db: 	[0.82,0,0.31],
				Sg: 	[0.85,0,0.27],
				Bh: 	[0.88,0,0.22],
				Hs: 	[0.9,0,0.18],
				Mt: 	[0.92,0,0.15]
			}, {
				get: function( target, name ){
						return target.hasOwnProperty( name ) ? target[name]: [1, 0.1, 0.55];
				}
			});

			this.bondCols = new Proxy( {
				1 : [0.1, 0.1, 0.1],
				2 : [0.5, 0.1, 0.3],
				3 : [0.1, 0.2, 0.5],
				4 : [0.9, 0.5, 0.2],
				9 : [0.3, 0.3, 0.3]
			}, {
				get: function( target, name ){
						return target.hasOwnProperty( name ) ? target[name]: [0, 0, 0];
				}
			});

			this.groupCols = new Proxy( {
				Carboxyl : [0.3, 0.3, 0.3],
				"Carboxyl + Alkyl (Ester)" : [0.3, 0.3, 0.3],
				Amide : [0.56, 0.56, 1.0],
				"Acyl Chloride" : [1.0, 1.0, 0.0],
				Acetal : [1.0, 0.0, 0.0],
				Carbonyl : [1.0, 0.5, 0.5],
				Alkoxy : [1.0, 0.0, 0.0],
				Amino : [0.56, 0.56, 1.0],
				Nitro : [0.56, 0.56, 1.0],
				Hydroxyl : [1.0, 0.0, 0.0],
				Fluoro : [0.0, 1.0, 0.0],
				Chloro : [0.0, 1.0, 0.0],
				Bromo : [0.64, 0.18, 0.18],
				Iodo : [0.0, 1.0, 0.0],
				Nitrile : [0.56, 0.56, 1.0],
			}, {
				get: function( target, name ){
						return target.hasOwnProperty( name ) ? target[name]: [1, 0.1, 0.55];
				}
			});
		}

		this._onWindowResize = function(){

			const elBox = self._DOM.getBoundingClientRect();
			self._aspect = elBox.width / elBox.height;
			self._camPersp.aspect = self._aspect;

			self.camActive.updateProjectionMatrix();
			self.renderer.setSize( elBox.width, elBox.height );
			self.controls.screen = { left: elBox.left, top: elBox.top, width: elBox.width, height: elBox.height };
			self.effect.render( self.scene, self.camActive );

		}

		this._animate = function(){

			Object.values( self.frameFunctions ).forEach( el => el.enabled && el.fn( self ) )

			self.effect.render( self.scene, self.camActive );

		}

	}

	/////3D canvas methods/////
	Object.assign( Mol3D.prototype, {

		init: function(){

			if( this._DOM ){

				//////CAMERAS//////
				const elBox = this._DOM.getBoundingClientRect();
				this._aspect = elBox.width / elBox.height;

				this._camPersp = new THREE.PerspectiveCamera( this._FOV, this._aspect, 1, 100);
				this.camActive = this._camPersp;
				this.controls = new THREE.OrbitControls( this._camPersp, this._DOM );
				this.controls.enablePan = false;
				this._camOrtho = new THREE.OrthographicCamera( 0, 0, 0, 0, -100, 100 ); //frustum applied later
				this._disableInteractions && this.controls.dispose()

				//////LIGHTS//////
				const directionalLight = new THREE.DirectionalLight( 0xffffff, 1 );
				directionalLight.position.set( 3, 3, 3 );
				directionalLight.name = "KeyLight";
				this.scene.add( directionalLight );

				//////RENDERER SETUP//////
				this.renderer = new THREE.WebGLRenderer( {antialias: true, alpha: true} );
				this.renderer.setPixelRatio( window.devicePixelRatio );
				this.renderer.setClearColor( 0xffffff, 0 );
				this.renderer.setSize( this._DOM.getBoundingClientRect().width, this._DOM.getBoundingClientRect().height );
				this._DOM.appendChild( this.renderer.domElement );

				this.effect = new THREE.OutlineEffect( this.renderer, {defaultThickness: 0.002} );

				this._initialised = true;

			}else{

				console.warn( "No container element specified! Set property '.Container = element'")

			}

			return this

		},

		draw: function() {

			if( !this._molecule ){ console.warn( "No molecule to draw! Set .Molecule = Molecule" ) }
			else if( this._molecule.ajaxRunning ){ console.warn( "AJAX Request currently in progress: no data to draw. Listen for event 'ajaxComplete' after calling getFromSMILE() to call draw()." ) }
			else if( !this._initialised ){ console.warn( "3D Container needs initialising. Call .init() before attempting to draw." ) }
			else{

				if (Detector.webgl){

					this.scene.remove( this.molGroup );
					this.molGroup = this.genGroup( this._molecule );
					this.scene.add( this.molGroup );
					this.showHs = this._showHs;
					this.resetView();

					if( this.animID ){ window.cancelAnimationFrame( this.animID ) };
					this.play();
					this._onWindowResize();
					window.addEventListener( "resize", this._onWindowResize );
					["mousemove","touchmove","touchstart"].forEach( e => document.addEventListener( e, evt => {

							//evt.preventDefault();

							const clientX = evt.changedTouches ? evt.changedTouches[0].clientX : evt.clientX;
							const clientY = evt.changedTouches ? evt.changedTouches[0].clientY : evt.clientY;


							const elBox = this.renderer.domElement.getBoundingClientRect()
							this.mouse.x = ( clientX - elBox.left ) / elBox.width*2 - 1;
							this.mouse.y = 1 - ( ( clientY - elBox.top ) / elBox.height*2);

						}, false )
					);

				}else{

					d3.select( this._DOM ).appendChild( Detector.getWebGLErrorMessage() );

				}
			}

			return this

		},

		dispose: function() {

			if( this.animID ){ this.pause() };

			this.scene.children.forEach( obj => this.scene.remove( obj ) );
			this.scene = null;

			this.renderer.domElement.remove();
			this.renderer.dispose();
			this.controls.dispose();

		},

		genGroup: function( molecule ){

			const group = new THREE.Group();

			group.add( ...this.drawAtoms( molecule.atoms ) );
			group.add( ...this.drawBonds( molecule.bonds ) );
			this.drawFunctionalGroups( molecule.fGroups, group );

			return group

		},

		drawAtoms: function( atoms ){

			var group = [];

			atoms.forEach( ( el, i ) => {

				if( !this.atomCols[el.element]["material"] ){

					this.atomCols[el.element]["material"] = new THREE.MeshToonMaterial({
							color: new THREE.Color().setRGB( ...this.atomCols[el.element] ),
							reflectivity: 0.8,
							shininess: 0.8,
							specular: 0.8,
						})

				}

				const sphere = new THREE.SphereGeometry( ( el.element === "H" ? 0.2 : 0.25 ) , 16, 16 )
				const mesh = new THREE.Mesh( sphere, this.atomCols[el.element]["material"] );
				mesh.position.set( ...el.pos );

				mesh.name = el.index;
				mesh.userData.tooltip = el.element;
				mesh.userData.type = "atom";
				mesh.userData.source = el;

				el.HTML = mesh;

				group.push( mesh );

			})

			return group
		},

		drawBonds: function( bonds ){

			var group = [];

			bonds.forEach( ( el, i ) => {

				var vec = new THREE.Vector3( ...el.end.pos ).sub( new THREE.Vector3( ...el.start.pos ) );

				switch ( el.type ){

					case 1:

						var bond = new THREE.TubeGeometry(
							new THREE.LineCurve3(
								new THREE.Vector3( 0, 0, 0 ),
								vec,
							), 0, .05, 14, false );
						break;

					case 2:
					case 3:

						var polar_sphere = [vec.length(), Math.acos(vec.z/vec.length()), Math.atan2(vec.y,vec.x)]; //[r,theta,phi]
						var bond = new THREE.Geometry();

						[0.5,-0.5].forEach( el => { //[Bend factor+, Bend factor-]

							bond.merge( new THREE.TubeGeometry(
								new THREE.QuadraticBezierCurve3(
									new THREE.Vector3( 0, 0, 0 ),
									new THREE.Vector3(
										vec.x/2 + el * Math.sin( polar_sphere[1] ) * Math.cos( polar_sphere[2] + Math.PI/2 ),
										vec.y/2 + el * Math.sin( polar_sphere[1] ) * Math.sin( polar_sphere[2] + Math.PI/2 ),
										vec.z/2 + el * Math.cos( polar_sphere[1] )
 									),
									vec
								), 10, .05, 8, false )
							)

						})

						if ( el.type === 3 ) {

							bond.merge(
								new THREE.TubeGeometry(
									new THREE.LineCurve3(
										new THREE.Vector3( 0, 0, 0 ),
										vec,
								), 0, .05, 14, false )
							);

						};
						break;

					case 4:

						var polar_sphere = [vec.length(), Math.acos(vec.z/vec.length()), Math.atan2(vec.y,vec.x)]; //[r,theta,phi]
						var bond = new THREE.Geometry();

						const offset = new THREE.Vector3(
							Math.sin( polar_sphere[1] )*Math.cos( polar_sphere[2] + Math.PI/2 ),
							Math.sin( polar_sphere[1] )*Math.sin( polar_sphere[2] + Math.PI/2 ),
							Math.cos( polar_sphere[1] )
						).multiplyScalar( 0.1 );

						bond.merge( new THREE.TubeGeometry(
							new THREE.LineCurve3(
								offset,
								new THREE.Vector3().copy( offset ).add( vec )
							), 0, .05, 8, false )
						);

						offset.negate();

						[[0.0, 0.15],[0.25, 0.45],[0.55, 0.75],[0.85, 1.0]].forEach( el => {
							bond.merge( new THREE.TubeGeometry(
								new THREE.LineCurve3(
									new THREE.Vector3().copy( offset ).addScaledVector( vec, el[0] ),
									new THREE.Vector3().copy( offset ).addScaledVector( vec, el[1] ),
									), 0, .05, 14, false )
							);
						});
						break;

					case 9:

						var bond = new THREE.Geometry();

						[[0.0, 0.15],[0.25, 0.45],[0.55, 0.75],[0.85, 1.0]].forEach( el => {
							bond.merge( new THREE.TubeGeometry(
								new THREE.LineCurve3(
									new THREE.Vector3().copy( vec ).multiplyScalar( el[0] ),
									new THREE.Vector3().copy( vec ).multiplyScalar( el[1] ),
									), 0, .05, 14, false )
							);
						})
						break;

				};

				if( !this.bondCols[el.type]["material"] ){

					this.bondCols[el.type]["material"] = new THREE.MeshToonMaterial({
							color: new THREE.Color().setRGB( ...this.bondCols[el.type] ),
							reflectivity: 0.8,
							shininess: 0.8,
							specular: 0.8,
						})

				}


				const mesh = new THREE.Mesh( bond, this.bondCols[el.type]["material"] );

				mesh.name = el.start.index + "_" + el.end.index;
				mesh.userData.source = el;
				mesh.userData.tooltip = el;
				mesh.userData.type = "bond";

				mesh.position.set( ...el.start.pos );

				el.HTML = mesh;

				group.push( mesh );

			})

			return group

		},

		drawFunctionalGroups: function( fGroups, parentGroup ){

			fGroups.forEach( function( el, i ){

				var hlmat = new THREE.MeshToonMaterial( {
							color: new THREE.Color().setHex( Math.random() * 0xffffff ),
							transparent: true,
							visible: true,
							opacity: 0.7,
							flatShading: true,
						});

				var hlmesh = new THREE.Mesh( new THREE.IcosahedronBufferGeometry( 0.5, 1 ), hlmat );

				hlmesh.userData.tooltip = el.type;
				hlmesh.userData.type = "fGroup";

				el.domain.forEach( function( fGroupDomain ){

					const hlmeshNEW = hlmesh.clone();
					hlmeshNEW.userData.source = el;
					parentGroup.getObjectByName( fGroupDomain.index ).add( hlmeshNEW );

					hlmeshNEW.visible = this._showfGroups

				});

				hlmesh.userData.source = el;

				parentGroup.getObjectByName( el.source.index ).add( hlmesh )

				hlmesh.visible = this._showfGroups

			})

		},

		resetView: function(){

			const boundSph = new THREE.Box3().setFromObject( this.molGroup ).getBoundingSphere();
			const angularSize = this.camActive.fov * Math.PI / 180;
			const distance = Math.sqrt( 2 ) * boundSph.radius / Math.tan( angularSize / 2 );

			this.setView( boundSph, new THREE.Vector3().set( distance, 0, 0), new THREE.Vector3( 0, 1, 0 ) );

			this.disableInteractions = this._disableInteractions;

		},

		printHierarchy: function( obj ) {

			if( !obj ){ obj = this.molGroup }

			function print( obj ){

				console.group( ' <' + obj.type + '> ' + obj.name );
				obj.children.forEach( self.printHierarchy );
				console.groupEnd();

			}

		},

		setView: function( lookAt, position, up ){

			this.camActive.lookAt( lookAt.center );
			this.camActive.position.copy( position );
			this.camActive.updateProjectionMatrix();
			this.camActive.up.copy( up )

			//this.controls.dispose()
			//this.controls = new THREE.OrbitControls( this.camActive, this._DOM );
			this.controls.update();

		},

		toggleCam: function(){

			this.camActive = this.camActive === this._camPersp ? this._camOrtho : this._camPersp;

		},

		pause: function(){

			cancelAnimationFrame( this.animID );
			this.animID = null;

		},

		play: function(){

			if( !this.animID ) this.animID = requestAnimationFrame( this._animate );

		},

	})

	/////3D canvas properties/////
	Object.defineProperties( Mol3D.prototype, {

		"disableInteractions": {

			get: function(){

				return this._disableInteractions;

			},

			set: function( value ){

				if( this.controls ){

					if( !value ){

						this.controls.enabled = true;
						d3.select( this._DOM ).style( "cursor", "all-scroll" );

					} else{

						this.controls.enabled = false;
						d3.select( this._DOM ).style( "cursor", null );

					}

				}

				this._disableInteractions = value;

			}

		},

		"showStats": {

			get: function(){

				return this._showStats;

			},

			set: function( value ){

				if( value ){

					this.stats = new Stats();
					this._DOM.appendChild( this.stats.dom )

				}else{

					if( this.stats ){ this.stats.dom.remove() };

				}

				this._showStats = value;

			}

		},

		"Container": {

			get: function(){

				return this._DOM;

			},

			set: function( value ){

				this._DOM = value;

			}

		},

		"Molecule": {

			get: function(){

				return this._molecule;

			},

			set: function( value ){

				this._molecule = value;

			}

		},

		"showHs": {

			get: function(){ return this._showHs },

			set: function( value ){

				this.molGroup.traverse(  ob => {

					switch( ob.userData.type ){

						case "atom":

							if( ob.userData.source.element === "H" && ob.userData.source.bondedTo[0].el.element === "C" ) ob.visible = value;
							break;

						case "bond":

							if( ( ob.userData.source.start.element === "H" && ob.userData.source.end.element === "C" ) || ( ob.userData.source.start.element === "C" && ob.userData.source.end.element === "H" ) ) ob.visible = value;
							break;

						case "fGroup":

							break;

					}
				})

			}

		},

		"showfGroups": {

			get: function(){ return this._showfGroups },

			set: function( value ){

				this.scene.traverse( obj => {
					if( obj.userData.type === "fGroup" ){
						obj.visible = value
					}
				})

				this._showfGroups = value;

			}

		},

		"onMouseDown": {

			get: function(){

				return this._onMouseDown

			},

			set: function( value ) {

				d3.select( this._DOM ).on( "mousedown", value, true );
				this._onMouseDown = value;

			}

		},

		"onWindowResize": {

			get: function(){

				return this._onWindowResize

			},

			set: function( value ){

					if( this._initialised ){ window.addEventListener( "resize", this._onWindowResize ) };
					this._onWindowResize = value;

			}

		},

		"highlight": {

			get: function(){ return this.frameFunctions["highlight"].enabled },

			set: function( value ){ this.frameFunctions["highlight"].enabled = value }

		},

		"autoRotate": {

			get: function(){ return this.frameFunctions["autoRotate"].enabled },

			set: function( value ){ this.frameFunctions["autoRotate"].enabled = value }

		},

		"highlightSync": {

			get: function(){ return this.frameFunctions["highlightSync"].enabled },

			set: function( value ){ this.frameFunctions["highlightSync"].enabled = value }

		},

		"mouseoverDispatch": {

			get: function(){ return this.frameFunctions["mouseoverDispatch"].enabled },

			set: function( value ){ this.frameFunctions["mouseoverDispatch"].enabled = value }

		},

		"labelTrack": {

			get: function(){ return this.frameFunctions["labelTrack"].enabled },

			set: function( value ){ this.frameFunctions["labelTrack"].enabled = value }

		},

		"labels": {

			get: function(){ return this.frameFunctions["labelTrack"].props.labels },

			set: function( value ){ this.frameFunctions["labelTrack"].props.labels = value }

		},

		"_frustum": {

			get: function(){

				let depth = this.controls.target.clone().sub( this._camPersp.position ).dot( this._camPersp.getWorldDirection( new THREE.Vector3() ) );
				let height_ortho = depth * 2 * Math.atan( this._FOV*( Math.PI/180 ) / 2 );
				let width_ortho = height_ortho * this._camPersp.aspect;

				return [width_ortho, height_ortho]

			},

			set: function( value ){}

		}

	})

	exports.Mol2D = Mol2D;
	exports.Mol3D = Mol3D;
	exports.Molecule = Molecule;
	exports.Atom = Atom;
	exports.Bond = Bond;

}))
