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
		this.object2D = null;
		this.object3D = null;

	}

	/////Bond object/////
	function Bond( index, bondStart, bondEnd, bondType, bondDirection ){

		this.index = index;
		this.start = bondStart;
		this.end = bondEnd;
		this.type = bondType;
		this.direction = bondDirection;
		this.claimed = false;
		this.object2D = null;
		this.object3D = null;

	}

	//////Functional Group object//////
	function fGroup( domain, claimed, type ){

		this.type = type;
		this.domain = domain;
		this.claimed = claimed;
		this.material = null;
		this.object3D = null;

	}

	/////Molecule object/////
	function Molecule( molFile ){

		if( molFile ){

			this.molFile = molFile;
			this.parseMol();

		} else{

			this.atoms = [];
			this.bonds = [];
			this.fGroups = [];
			this.molFile = null;

		}

	};

	/////Molecule methods/////
	Object.assign( Molecule.prototype, {

		parseMol: function(){

			let mol = {};
			molFile = this.molFile.split( "\n" ).map( el => el.trim() ).join( "\n")
			const [numatoms,numbonds] = molFile.split( "\n" )[3].match( /.{1,3}/g ).slice( 0, 2 ).map( el => parseInt( el ) );

			mol.atoms = molFile.split( "\n" )
							.slice( 4, 4 + numatoms )
							.map( function( el, i ){

								el = " ".repeat( 69 - el.length ) + el;
								const line = el.slice( 0, 30 ).match( /.{1,10}/g ).concat( el.slice( 30 ).match( /.{1,3}/g ) );
								return new Atom(
										i,
										{ x: parseFloat( line[0] ), y: parseFloat( line[1] ), z: parseFloat( line[2] ) },
										line[3].trim(),
										0
									)

							});
			mol.bonds = molFile.split( "\n" )
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
			molFile.split( "\n" ).forEach( function( el, i ){
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

			const fGroupSearcher = () => {

				let fGroups = []

				mol.atoms.filter( atom => atom.element != "H" ).forEach( atom => {
					let scanning = true

					while( scanning ){
						scanning = false

						for( ss of MolViewer.subStructure_Lib.filter( ss => ( ss.root === "R" ? true : ss.root === atom.element ) ) ) { //Filter out non-matching root atoms

							let foundStructures = inBonds( atom, ss.bonds, atom.index, [] );

							if( foundStructures ){
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
					let results;

					subStruct.forEach( function( ss ){

						ss.found = false;
						source.bondedTo
							.filter( bond => bond.el.index !== rootIndex && !domain.map( dom => dom.index ).includes( bond.el.index ) && bond.bond.claimed === false )
							.forEach( function( bond ){

								if( !ss.found && bond.bond.type === ss.btype ){
									if( ( ss.el === "R" ? true : ( ss.el === "X" ? ["Cl", "Br", "I", "F"].includes( bond.el.element ) : bond.el.element === ss.el ) ) ){

										domain.push( bond.el );
										claimed.push( bond.bond );

										if( ss.hasOwnProperty( "bondedTo" ) ){

											//////Recursive search//////
											const deepSearch = inBonds( bond.el , ss.bondedTo, rootIndex, domain );

											if( deepSearch ){
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
						results = new MolViewer.fGroup( [source].concat( domain ), claimed );
					} else{
						results = null;
					}

					return results

				}

				return fGroups

			}

			this.atoms = mol.atoms;
			this.bonds = mol.bonds;
			this.fGroups = fGroupSearcher();
			this._molFile = molFile;

			return this;
		},

		get2DFromSMILE: function( smile, addHydrogens ){

			let molecule = OCL.Molecule.fromSmiles( smile );
			addHydrogens && molecule.addImplicitHydrogens();
			this.molFile = molecule.toMolfile();
			this.parseMol();

			return this

		},

		get3DFromSMILE: function( smile ){

			this.ajaxRunning = true;
			const event = new Event( "ajaxComplete" );

			d3.xml( "https://cactus.nci.nih.gov/chemical/structure/" + smile.replace( /\#/g, "%23" ).replace( /\[/g, "%5B" ).replace( /\]/g, "%5D" ) + "/file/xml?format=sdf&get3d=true")
				.then( ( d, err ) => {

					if( d.children[0].children.length === 0 ){ console.error( "MolViewer: Invalid SMILE string!"); return };

					this.molFile = d.children[0].children[0].children[0].innerHTML;
					this.parseMol();
					this.ajaxRunning = false
					document.dispatchEvent( event )

				} )

		},

		centre: function(){

			const centre = this.atoms.reduce( ( acc, atom ) => { acc[0] += atom.pos.x; acc[1] += atom.pos.y; acc[2] += atom.pos.z; return acc }, [0,0,0] ).map( coord => coord/this.atoms.length );

			this.atoms = this.atoms.map( atom => { atom.pos.x -= centre[0]; atom.pos.y -= centre[1]; atom.pos.z -= centre[2]; return atom } );

			return this;

		}

	})

	/////Molecule properties/////
	Object.defineProperties( Molecule.prototype , {

		"bondLength": {

		    get: function(){

				return this.bonds.reduce( ( acc, bond ) =>
					acc + Math.sqrt( Math.pow( bond.start.pos.x - bond.end.pos.x, 2 ) +
									 Math.pow( bond.start.pos.y - bond.end.pos.y, 2 ) +
									 Math.pow( bond.start.pos.z - bond.end.pos.z, 2 ) )
					, 0 ) / this.bonds.length

			},

		    set: function( value ){

				const scale = value/this.bondLength;

				this.atoms = this.atoms.map( atom => { atom.pos.x *= scale; atom.pos.y *= scale; atom.pos.z *= scale; return atom } );

			}

		},

	})

	/////2D canvas objects/////
	function Mol2D( Molecule, Container, dims, params ){

		if( Container === undefined ){ console.error( "No container element specified!" ); return };

		params = params || {};

		const self = this;

		this._initialised = false;
		this._bondScale = 50;
		this._scaleBox = { width : 0, height: 0 };
		this.dims = dims !== undefined ? dims : { x: 0, y: 0, width: Container.getBoundingClientRect().width, height: Container.getBoundingClientRect().height };
		this.Container = Container;
		this.Molecule = Molecule;

		this.zoomable = params.zoomable !== undefined ? params.zoomable : true;
		this.showIndices = params.showIndices !== undefined ? params.showIndices : false;
		this.showHs = params.showHs !== undefined ? params.showHs : true;
		this.zoomSmooth = params.zoomSmooth !== undefined ? params.zoomSmooth : 300;
		this.zoomEase = params.zoomEase !== undefined ? params.zoomEase : d3.easeCircleOut;

		this.stylesheet = `
			#view2d {
				text-anchor: middle;
				font-family: sans-serif;
				font-size: 16px;
				cursor: all-scroll;
			}

			.bondline, .bond_dbl > line, .bond_hash > line, .bond_trp > line {
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
		`;

		this._onWindowResize = function( ev ){

			if( self.svg ){

				const contBox = self._DOM.getBoundingClientRect();

				self.dims.width = self._scaleBox.width * contBox.width;
				self.dims.height = self._scaleBox.height * contBox.height;

			}

		}

	};

	/////2D canvas methods/////
	Object.assign( Mol2D.prototype, {

		init: function(){

			if( this._initialised ) console.warn( "MolViewer: Molecule already initialised!" );

			d3.select( "#molViewer2DCSS" ).remove();
			d3.select( "html" ).append( "style" ).attr( "id", "molViewer2DCSS" ).html( this.stylesheet );
			this._initialised = true;

			return this

		},

		draw: function(){

			if( this.svg ){ this.svg.remove(); this.svg = null };
			const svg = this.Container.append( "svg" ).attr( "viewBox", this.dims.x + " " + this.dims.y + " " + this.dims.width + " " + this.dims.height ).attr( "id", "view2d" );
			this.dims = svg.node().viewBox.baseVal;

			this.svg = svg;
			this.root = this.genMolecule();
			this.svg.node().appendChild( this.root.node() );
			this.showHs = this.showHs;

			this.zoomable && this.fitToScreen();

			window.removeEventListener( "resize", this._onWindowResize );
			window.addEventListener( "resize", this._onWindowResize );

			return this

		},

		genMolecule: function(){

			const atoms = this.Molecule.atoms;
			const bonds = this.Molecule.bonds;
			const bondScale = this.bondScale;

			const root = d3.create( "svg:g" )
					.attr( "id", "rootframe" )
					.attr( "transform", null )

			//////ATOMS//////
			root.append( "g" )
				.attr( "class", "atoms" )
				.selectAll( "g" )
				.data( this.Molecule.atoms )
				.enter()
				.each( drawAtom )

			//////BONDS//////
			const labelOffset = 8;
			const bondsroot = root.append( "g" )
				.attr( "class", "bonds" )
				.selectAll( "g" )
				.data( this.Molecule.bonds )
				.enter()
				.each( drawBond )

			function drawAtom( atom, i ){

				const pos = { x: atom.pos.x * bondScale, y: atom.pos.y * -bondScale };

				const atomGrp = d3.select( this ).append( "g" )
									.attr( "class", "atom_" + atom.element )
									.attr( "id", i )

				atom.element != "H" && atomGrp.append( "g" ).attr( "class", "hydrogens" );

				const highlightCircle = atomGrp.append( "circle" )
									.attr( "class", "highlight" )
									.attr( "id", "highlight_" + i )
									.attr( "cx", pos.x )
									.attr( "cy", pos.y )
									.attr( "r", bondScale/3 )

				const txt = atomGrp.append( "text" )
									.attr( "class", "label_" + atom.element )
									.attr( "id", "label_" + i )
									.attr( "x", pos.x )
									.attr( "y", pos.y + 6 )
									.text( atom.element !== "C" || atom.charge ? atom.element : "" )

				const ind = atomGrp.append( "g" )
									.attr( "class", "atomind" )
									.attr( "id", "atomind_" + i )
									.attr( "display", this.showIndices ? null : "none" )

				const indBBox = ind.append( "text" )
									.attr( "x", pos.x - txt.node().getBBox().width/2 - i.toString().split("").length * 2 )
									.attr( "y", pos.y - 5 )
									.text( i )
									.node().getBBox()

				ind.append( "rect" )
					.attr( "x", indBBox.x )
					.attr( "width", indBBox.width )
					.attr( "y", indBBox.y )
					.attr( "height", indBBox.height )
					.lower();

				if( atom.charge && false ){ //TODO

					const chg = atomGrp.append( "g" )
									.attr( "class", "charge" );
					const chgTxt = chg.append( "text" )
										.attr( "class", "atomcharge" )
										.attr( "id", "atomchg_" + i )
										.attr( "x", atom.pos.x + txt.node().getBBox().width/2 + atom.charge.toString().length * 2 - 2 )
										.attr( "y", atom.pos.y - 8 )
										.text( atom.charge === -1 ? "-" : ( atom.charge === 1 ? "+" : atom.charge ) )

					chg.append( "circle" )
						.attr( "cx", chgTxt.attr( "x" ) )
						.attr( "cy", chgTxt.attr( "y") )
						.attr( "r", chgTxt.node().getBBox().width/2 + 1 )
						.lower();

				}

				atom.object2D = atomGrp.node();

			}

			function drawBond( bond, i ){

				let tmp;

				const posStart = { x: bond.start.pos.x * bondScale, y: bond.start.pos.y * -bondScale };
				const posEnd = { x: bond.end.pos.x * bondScale, y: bond.end.pos.y * -bondScale };

				const bondGrp = d3.select( this ).append( "g" ).attr( "class", "bond_" + bond.type );
				const theta = Math.atan2( posEnd.y - posStart.y, posEnd.x - posStart.x );

				const coords = [posStart.x + ( bond.start.element !== "C" || bond.start.charge ? bond.start.element.length * labelOffset * Math.cos( theta ) : 0 ),
								posEnd.x - ( bond.end.element !== "C" || bond.end.charge ? bond.end.element.length * labelOffset * Math.cos( theta ) : 0 ),
								posStart.y + ( bond.start.element !== "C" || bond.start.charge ? bond.start.element.length * labelOffset * Math.sin( theta ) : 0 ),
								posEnd.y - ( bond.end.element !== "C" || bond.end.charge ? bond.end.element.length * labelOffset * Math.sin( theta ) : 0 )];

				const length = Math.hypot( coords[1] - coords[0], coords[2] - coords[3] );

				const highlight = bondGrp.append( "rect" )
					.attr( "x", -bondScale/6 ).attr( "y", -bondScale/6 )
					.attr( "rx", bondScale/10 ).attr( "ry", bondScale/10 )
					.attr( "width", length + bondScale/3 ).attr( "height", bondScale/3)
					.attr("transform", "translate(" + posStart.x + "," + posStart.y + ")rotate(" + theta*180/Math.PI + ")" )
					.attr( "class", "highlight" )
					.attr( "id", "highlight_" + bond.start.index + "_" + bond.end.index );


				const placeholderLine = bondGrp.append( "line" )
						.attr( "class", "bondline")
						.attr( "x1", coords[0] )
						.attr( "x2", coords[1] )
						.attr( "y1", coords[2] )
						.attr( "y2", coords[3] )

				switch( bond.type ){

					case 1: //single bond

						switch( bond.direction ){

							case 0: //normal bond

								bondGrp.attr( "class", "bond" );
								break;

							case 1: //wedge bond

								tmp = bondGrp.attr( "class", "bond_wedge" ).append( "polygon" )
									.attr( "points", coords[0] + "," + coords[2] + " " +
										 ( coords[1] + 3*Math.cos( theta + Math.PI/2 ) ).toFixed( 2 ) + "," + ( coords[3] + 3*Math.sin( theta + Math.PI/2 ) ).toFixed( 2 ) + " " +
										 ( coords[1] - 3*Math.cos( theta + Math.PI/2 ) ).toFixed( 2 ) + "," + ( coords[3] - 3*Math.sin( theta + Math.PI/2 ) ).toFixed( 2 ) );
								placeholderLine.remove()
								break;

							case 6: //hash bond

								tmp = bondGrp.attr( "class", "bond_hash" );
								let point = [0, 0];
								for( let i = 0; Math.hypot( ...point ) < length ; i++ ){
									tmp.append( "line" ).attr( "class", "bond" )
										.attr( "x1", ( coords[0] + point[0] + Math.hypot( ...point )/length * 3*Math.cos( theta + Math.PI/2 ) ).toFixed( 2 ) )
										.attr( "x2", ( coords[0] + point[0] - Math.hypot( ...point )/length * 3*Math.cos( theta + Math.PI/2 ) ).toFixed( 2 ) )
										.attr( "y1", ( coords[2] + point[1] + Math.hypot( ...point )/length * 3*Math.sin( theta + Math.PI/2 ) ).toFixed( 2 ) )
										.attr( "y2", ( coords[2] + point[1] - Math.hypot( ...point )/length * 3*Math.sin( theta + Math.PI/2 ) ).toFixed( 2 ) );
									point = [( i + 1 ) * 3 * Math.cos( theta ), ( i + 1 ) * 3 * Math.sin( theta )];
								};
								placeholderLine.remove()
								break;

						}

						break;

					case 2: //double bond

						bondGrp.attr( "class", "bond_dbl" );
						[-1, 1].forEach( el => {

							bondGrp.node().appendChild( d3.select( placeholderLine.node().cloneNode() )
								.attr( "x1", coords[0] + el*bondScale/15*Math.cos( theta + Math.PI/2 ) )
								.attr( "x2", coords[1] + el*bondScale/15*Math.cos( theta + Math.PI/2 ) )
								.attr( "y1", coords[2] + el*bondScale/15*Math.sin( theta + Math.PI/2 ) )
								.attr( "y2", coords[3] + el*bondScale/15*Math.sin( theta + Math.PI/2 ) )
								.node() );

						})
						placeholderLine.remove()
						break;

					case 3: //triple bond

						bondGrp.attr( "class", "bond_trp" );
						[-1, 1].forEach( el => {

							bondGrp.node().appendChild( d3.select( placeholderLine.node().cloneNode() )
								.attr( "x1", coords[0] + el*bondScale/15*Math.cos( theta + Math.PI/2 ) )
								.attr( "x2", coords[1] + el*bondScale/15*Math.cos( theta + Math.PI/2 ) )
								.attr( "y1", coords[2] + el*bondScale/15*Math.sin( theta + Math.PI/2 ) )
								.attr( "y2", coords[3] + el*bondScale/15*Math.sin( theta + Math.PI/2 ) )
								.node() );

						})
						break;

					case 9: //aromatic bond

						placeholderLine.attr("class","bond")
						break;

				}

				if( bond.start.element === "H" || bond.end.element === "H" ){

					const rootAtom = root.select( "[id='" + ( bond.start.element === "H" ? bond.end.index : bond.start.index ) + "'] > .hydrogens" ).node();

					rootAtom.appendChild( bondGrp.node() );
					rootAtom.appendChild( root.select( "[id='" + ( bond.start.element !== "H" ? bond.end.index : bond.start.index ) + "']" ).node() );

				}

				bond.object2D = bondGrp.node();

			}

			return root

		},

		fitToScreen: function(){

			const zoomfn = d3.zoom().on( "zoom", () => {

				this.root.transition()
					.duration( this.zoomSmooth )
					.ease( this.zoomEase )
					.attr( "transform", d3.event.transform )

			})

			this.root.attr( "transform", null ) ;
			const viewBox = this.root.node().parentNode.viewBox.baseVal;
			const rootBox = this.root.node().getBBox();

			const zoom = viewBox.width/rootBox.width < viewBox.height/rootBox.height ? viewBox.width/rootBox.width : viewBox.height/rootBox.height;

			this.svg.call( zoomfn.transform, d3.zoomIdentity.translate( viewBox.width/2 + ( - ( rootBox.x - viewBox.x ) - rootBox.width/2 )*zoom, viewBox.height/2 + ( - ( rootBox.y - viewBox.y ) - rootBox.height/2 )*zoom ).scale( zoom ) )
			this.svg.call( zoomfn );

		},

	})

	/////2D canvas properties/////
	Object.defineProperties( Mol2D.prototype, {

		"Container": {

			get: function(){

				return d3.select( this._DOM );

			},

			set: function( value ){

				if( value instanceof Element ){

					this._DOM = value;
					this._scaleBox = { width: this.dims.width/value.getBoundingClientRect().width, height: this.dims.height/value.getBoundingClientRect().height };

				} else{

						console.warn( "MolViewer: Container is not of type Element!" );

				}

			}

		},

		"zoomable": {

			get: function () {

				return this._zoomable;

			},

			set: function ( value ) {

				this._zoomable = value;

			}

		},

		"zoomSmooth": {

		    get: function(){ return this._zoomSmooth },

		    set: function( value ){

				if( +value === +value ){

					this._zoomSmooth = value;

				} else{

					console.warn( "MolViewer: zoomSmooth must be an integer!" );

				}

			}

		},

		"zoomEase": {

		    get: function(){ return this._zoomEase },

		    set: function( value ){ this._zoomEase = value }

		},

		"showIndices": {

			get: function () {

				return this._showIndices;

			},

			set: function ( value ) {

				this._showIndices = value;
				if( !this.root ) return;

				this.root.selectAll( ".atomind" ).attr( "display", value ? null : "none" );

			}

		},

		"bondScale": {

		    get: function(){ return this._bondScale },

		    set: function( value ){ this._bondScale = value; this.draw(); }

		},

		"showHs": {

			get: function(){ return this._showHs },

			set: function( value ){

				this._showHs = value;
				if( !this.root ) return;

				this.root.selectAll( ".hydrogens, .hydrogens > *" ).each( function(){

					d3.select( this ).attr( "display", value ? null : "none" );

				})


			}

		},

	})

	/////3D canvas object/////
	function Mol3D( molecule, container, params ){

		params = params || {};

		const self = this;

		this.Container    = container || null;
		this.Molecule     = molecule;
		this.Scene        = new THREE.Scene();
		this.mouse        = new THREE.Vector2();
		this.molGroup     = new THREE.Group();
		this._animID;
		this._initialised = false;
		this._FOV = 70;

		this.stylesheet = `
			.view3D {
				text-anchor: middle;
				font-family: sans-serif;
				font-size: 16px;
			}
		`

		this.frameFunctions = {
			//////Scene controller (ONLY TOUCH IF YOU KNOW WHAT YOU'RE DOING)//////
			sceneController: {enabled: true, props: {}, fn: function( self ){

					self._animID = requestAnimationFrame( self._animate );

					self.effect.render( self.Scene, self.camActive );

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
			highlight: {enabled: false, props: { _intersected: null, _storedMat: null }, fn: function( self ){

					var raycaster = new THREE.Raycaster();
					raycaster.setFromCamera( self.mouse, self.camActive );

					const intersects = raycaster.intersectObjects( self.molGroup.children, true );

					if( self.mouse.x > -1 && self.mouse.x < 1 && self.mouse.y > -1 && self.mouse.y < 1 ){

						if( intersects.length > 0 ){

							if( intersects[0].object instanceof THREE.Mesh ){

								if( this.props._intersected !== intersects[0].object ){

									setMaterial( this.props._intersected, this.props._storedMat );

									this.props._storedMat = null;

								}

								this.props._intersected = intersects[0].object;
								if( !this.props._storedMat ) this.props._storedMat = this.props._intersected.material;
								const mat = this.props._storedMat.clone()
								mat.emissive.setHex( 0xff0000 )

								setMaterial( this.props._intersected, mat );

							}


						} else{

							setMaterial( this.props._intersected, this.props._storedMat );

							this.props._storedMat = null;

						}

					}

					function setMaterial( obj, mat ){

						if( !obj || !mat ) return;

						obj.material = mat;

						if( obj.userData.source instanceof MolViewer.fGroup ){

							obj.userData.source.object3D.forEach( atom => atom.material = mat );

						}

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

								self.hovered = d3.selectAll( highlighted.userData.source.domain.map( el =>  "#highlight_" + el.index ).join( ", " ) )
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
			mouseoverDispatch: {enabled: false, props: { _hovered: null }, fn: function( self ){

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

		this.highlight           = params.highlight !== undefined ? params.highlight : true;
		this.autoRotate          = params.autoRotate !== undefined ? params.autoRotate : false;
		this.highlightSync       = params.highlightSync !== undefined ? params.highlightSync : false;
		this.mouseoverDispatch   = params.mouseoverDispatch !== undefined ? params.mouseoverDispatch : false;
		this.labelTrack          = params.labelTrack !== undefined ? params.labelTrack : true;

		this.disableInteractions = params.disableInteractions !== undefined ? params.disableInteractions : false;
		this.showfGroups         = params.showfGroups !== undefined ? params.showfGroups : true;
		this.showHs              = params.showHs !== undefined ? params.showHs : true;
		this.showStats           = params.showStats !== undefined ? params.showStats : false;

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
						if( !target.hasOwnProperty( name ) ){ target[name] = [ 1, 0.1, 0.55 ] };
						return target[name]
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
					if( !target.hasOwnProperty( name ) ){ target[name] = [ 0, 0, 0 ] };
					return target[name]
				}
			});

			this.groupCols = new Proxy( {
				"Carboxylic Acid" : [0.3, 0.3, 0.3],
				Ester : [0.3, 0.3, 0.3],
				Amide : [0.56, 0.56, 1.0],
				"Acyl Halide" : [1.0, 1.0, 0.0],
				Aldehyde : [1.0, 0.0, 0.0],
				Ketone : [1.0, 0.5, 0.5],
				"Primary Amine" : [0.56, 0.56, 1.0],
				"Secondary Amine" : [0.36, 0.36, 0.8],
				"Tertiary Amine" : [0.16, 0.16, 0.6],
				Nitro : [0.56, 0.56, 1.0],
				Alcohol : [1.0, 0.0, 0.0],
				Halo : [0.0, 1.0, 0.0],
				Nitrile : [0.56, 0.56, 1.0],
				Acetal : [0.56, 0.56, 1.0],
				Ether : [0.56, 0.56, 1.0],
			}, {
				get: function( target, name ){
					if( !target.hasOwnProperty( name ) ){ target[name] = [1, 0.1, 0.55] };
					return target[name]
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
			self.effect.render( self.Scene, self.camActive );

		}

		this._animate = function(){

			Object.values( self.frameFunctions ).forEach( el => el.enabled && el.fn( self ) )

		}

	}

	/////3D canvas methods/////
	Object.assign( Mol3D.prototype, {

		init: function(){

			if( this._initialised ){ console.warn( "MolViewer: Molecule already initialised!" ); return }

			d3.select( "#molViewer3DCSS" ).remove();
			d3.select( "html" ).append( "style" ).attr( "id", "molViewer3DCSS" ).html( this.stylesheet );

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
				this.Scene.add( directionalLight );

				//////RENDERER SETUP//////
				this.renderer = new THREE.WebGLRenderer( {antialias: true, alpha: true} );
				this.renderer.setPixelRatio( window.devicePixelRatio );
				this.renderer.setClearColor( 0xffffff, 0 );
				this.renderer.setSize( this._DOM.getBoundingClientRect().width, this._DOM.getBoundingClientRect().height );
				this._DOM.appendChild( this.renderer.domElement );

				this.effect = new THREE.OutlineEffect( this.renderer, {defaultThickness: 0.002} );

				this._initialised = true;

			}else{

				console.warn( "MolViewer: No container element specified! Set property '.Container = element'")

			}

			return this

		},

		draw: function() {

			const updateMousePos = evt => {

				const clientX = evt.changedTouches ? evt.changedTouches[0].clientX : evt.clientX;
				const clientY = evt.changedTouches ? evt.changedTouches[0].clientY : evt.clientY;


				const elBox = this.renderer.domElement.getBoundingClientRect()
				this.mouse.x = ( clientX - elBox.left ) / elBox.width*2 - 1;
				this.mouse.y = 1 - ( ( clientY - elBox.top ) / elBox.height*2);

			}

			if( !this._molecule ){ console.warn( "MolViewer: No molecule to draw! Set .Molecule = Molecule" ) }
			else if( this._molecule.ajaxRunning ){ console.warn( "MolViewer: AJAX Request currently in progress: no data to draw. Listen for event 'ajaxComplete' after calling getFromSMILE() to call draw()." ) }
			else if( !this._initialised ){ console.warn( "MolViewer: 3D Container needs initialising. Call .init() before attempting to draw." ) }
			else{

				this.Scene.remove( this.molGroup );
				this.molGroup = this.genGroup( this._molecule );
				this.Scene.add( this.molGroup );
				this.showHs = this.showHs;
				this.resetView();

				if( this._animID ){ window.cancelAnimationFrame( this._animID ); this._animID = null; };
				this.play();
				this._onWindowResize();
				window.removeEventListener( "resize", this._onWindowResize );
				window.addEventListener( "resize", this._onWindowResize );
				["mousemove","touchmove","touchstart"].forEach( e => document.addEventListener( e, updateMousePos, false )
				);

			}

			return this

		},

		dispose: function() {

			if( this._animID ){ this.pause() };

			this.Scene.children.forEach( obj => this.Scene.remove( obj ) );
			this.Scene = null;

			this.renderer.domElement.remove();
			this.renderer.dispose();
			this.controls.dispose();

		},

		genGroup: function(){

			const drawAtom = ( parentGroup, atom, i ) => {

				if( !this.atomCols[atom.element]["material"] ){

					if( !this.atomCols[atom.element] ) this.atomCols[atom.element] = this.atomCols.xx; //Set it as default colour

					this.atomCols[atom.element]["material"] = new THREE.MeshToonMaterial({
							color: new THREE.Color().setRGB( ...this.atomCols[atom.element] ),
							reflectivity: 0.8,
							shininess: 0.8,
							specular: 0.8,
						})

				}

				const sphere = new THREE.SphereGeometry( ( atom.element === "H" ? 0.2 : 0.25 ) , 16, 16 )
				const mesh = new THREE.Mesh( sphere, this.atomCols[atom.element]["material"] );
				mesh.position.copy( atom.pos );

				mesh.name = atom.index;
				mesh.userData.tooltip = atom.element;
				mesh.userData.type = "atom";
				mesh.userData.source = atom;

				atom.object3D = mesh;

				parentGroup.add( mesh );

			}

			const drawBond = ( parentGroup, bond, i ) => {

				const vec = new THREE.Vector3().copy( bond.end.pos ).sub( new THREE.Vector3().copy( bond.start.pos ) );

				switch ( bond.type ){

					case 1:

						var bondGeo = new THREE.TubeGeometry(
							new THREE.LineCurve3(
								new THREE.Vector3( 0, 0, 0 ),
								vec,
							), 0, .05, 14, false );
						break;

					case 2:
					case 3:

						var polar_sphere = [vec.length(), Math.acos(vec.z/vec.length()), Math.atan2(vec.y,vec.x)]; //[r,theta,phi]
						var bondGeo = new THREE.Geometry();

						[0.5,-0.5].forEach( bend => { //[Bend factor+, Bend factor-]

							bondGeo.merge( new THREE.TubeGeometry(
								new THREE.QuadraticBezierCurve3(
									new THREE.Vector3( 0, 0, 0 ),
									new THREE.Vector3(
										vec.x/2 + bend * Math.sin( polar_sphere[1] ) * Math.cos( polar_sphere[2] + Math.PI/2 ),
										vec.y/2 + bend * Math.sin( polar_sphere[1] ) * Math.sin( polar_sphere[2] + Math.PI/2 ),
										vec.z/2 + bend * Math.cos( polar_sphere[1] )
 									),
									vec
								), 10, .05, 8, false )
							)

						})

						if ( bond.type === 3 ) {

							bondGeo.merge(
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
						var bondGeo = new THREE.Geometry();

						const offset = new THREE.Vector3(
							Math.sin( polar_sphere[1] )*Math.cos( polar_sphere[2] + Math.PI/2 ),
							Math.sin( polar_sphere[1] )*Math.sin( polar_sphere[2] + Math.PI/2 ),
							Math.cos( polar_sphere[1] )
						).multiplyScalar( 0.1 );

						bondGeo.merge( new THREE.TubeGeometry(
							new THREE.LineCurve3(
								offset,
								new THREE.Vector3().copy( offset ).add( vec )
							), 0, .05, 8, false )
						);

						offset.negate();

						[[0.0, 0.15],[0.25, 0.45],[0.55, 0.75],[0.85, 1.0]].forEach( tubeFragment => {
							bondGeo.merge( new THREE.TubeGeometry(
								new THREE.LineCurve3(
									new THREE.Vector3().copy( offset ).addScaledVector( vec, tubeFragment[0] ),
									new THREE.Vector3().copy( offset ).addScaledVector( vec, tubeFragment[1] ),
									), 0, .05, 14, false )
							);
						});
						break;

					case 9:

						var bondGeo = new THREE.Geometry();

						[[0.0, 0.15],[0.25, 0.45],[0.55, 0.75],[0.85, 1.0]].forEach( tubeFragment => {
							bondGeo.merge( new THREE.TubeGeometry(
								new THREE.LineCurve3(
									new THREE.Vector3().copy( vec ).multiplyScalar( tubeFragment[0] ),
									new THREE.Vector3().copy( vec ).multiplyScalar( tubeFragment[1] ),
									), 0, .05, 14, false )
							);
						})
						break;

				};

				if( !this.bondCols[bond.type]["material"] ){

					if( !this.bondCols[bond.type] ) this.bondCols[bond.type] = this.bondCols.xx; //Set it as default colour

					this.bondCols[bond.type]["material"] = new THREE.MeshToonMaterial({
							color: new THREE.Color().setRGB( ...this.bondCols[bond.type] ),
							reflectivity: 0.8,
							shininess: 0.8,
							specular: 0.8,
						})

				}

				const mesh = new THREE.Mesh( bondGeo, this.bondCols[bond.type]["material"] );

				mesh.name = bond.start.index + "_" + bond.end.index;
				mesh.userData.source = bond;
				mesh.userData.tooltip = bond;
				mesh.userData.type = "bond";

				mesh.position.copy( bond.start.pos );

				bond.object3D = mesh;

				parentGroup.add( mesh );

			}

			const drawFunctionalGroup = ( parentGroup, fGroup, i ) => {

				const fGroups = [];

				if( !fGroup.material ){

					if( !this.groupCols[fGroup.type]["material"] ){

						if( !this.groupCols[fGroup.type] ) this.groupCols[fGroup.type] = this.groupCols.xx; //Set as default colour

							this.groupCols[fGroup.type]["material"] = new THREE.MeshToonMaterial({
									color: new THREE.Color().setRGB( ...this.groupCols[fGroup.type] ),
									transparent: true,
									visible: true,
									opacity: 0.7,
									flatShading: true,
								})

					}

				}

				fGroup.material = this.groupCols[fGroup.type]["material"];

				const hlgeo = new THREE.IcosahedronBufferGeometry( 0.5, 1 );

				fGroup.domain.forEach( atom => {

					const hlmesh = new THREE.Mesh( hlgeo.clone(), fGroup.material );
					hlmesh.userData.source = fGroup;
					hlmesh.userData.tooltip = fGroup.type;
					hlmesh.userData.type = "fGroup";
					parentGroup.getObjectByName( atom.index ).add( hlmesh );
					fGroups.push( hlmesh );

					hlmesh.visible = this.showfGroups

				});

				fGroup.object3D = fGroups;

			};

			const group = new THREE.Group();

			this.Molecule.atoms.forEach( ( atom, i ) => drawAtom( group, atom, i ) );
			this.Molecule.bonds.forEach( ( bond, i ) => drawBond( group, bond, i ) );
			this.Molecule.fGroups.forEach( ( fGroup, i ) => drawFunctionalGroup( group, fGroup, i ) ); //Gets added directly to atoms, doesn't need adding to entire group

			return group

		},

		resetView: function(){

			const boundSph = new THREE.Box3().setFromObject( this.molGroup ).getBoundingSphere( new THREE.Sphere() );
			const angularSize = this.camActive.fov * Math.PI / 180;
			const distance = Math.sqrt( 2 ) * boundSph.radius / Math.tan( angularSize / 2 );

			this.setView( boundSph, new THREE.Vector3().set( 0, 0, distance), new THREE.Vector3( 0, 1, 0 ) );

			this.disableInteractions = this._disableInteractions;

		},

		printHierarchy: function( object ) {

			if( !object ){ object = this.molGroup }

			(function print( object ){

				console.group( ' <' + object.type + '> ' + object.name );
				object.children.forEach( print );
				console.groupEnd();

			})( object )

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

			cancelAnimationFrame( this._animID );
			this._animID = null;

		},

		play: function(){

			if( !this._animID ) this._animID = requestAnimationFrame( this._animate );

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

				if( value instanceof Element ){

					this._DOM = value;

				} else{

						console.warn( "MolViewer: Container is not of type Element!" );

				}

			}

		},

		"Molecule": {

			get: function(){

				return this._molecule;

			},

			set: function( value ){

				this._molecule = value instanceof MolViewer.Molecule ? value : new MolViewer.Molecule( value );

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

				this._showHs = value;

			}

		},

		"showfGroups": {

			get: function(){ return this._showfGroups },

			set: function( value ){

				this.Scene.traverse( obj => {
					if( obj.userData.type === "fGroup" ){
						obj.visible = value
					}
				})

				this._showfGroups = value;

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
	exports.subStructure_Lib = this.subStructures = [
						  {type: "Carboxylic Acid", root: "C", bonds:[{el: "O", btype: 2}, {el: "O", btype: 1, bondedTo:[{el: "H",  btype: 1}]}]},
						  {type: "Ester", root: "C", bonds:[{el: "O", btype: 2}, {el: "O", btype: 1, bondedTo:[{el: "R", btype:  1}]}]},
						  {type: "Amide", root: "C", bonds:[{el: "O", btype: 2}, {el: "N", btype: 1, bondedTo:[{el: "R", btype: 1},  {el: "R", btype: 1}]}]},

						  {type: "Acyl Halide", root: "C", bonds:[{el: "O", btype: 2}, {el: "X", btype: 1}]},

						  {type: "Aldehyde", root: "C", bonds:[{el: "O", btype: 2}, {el: "H", btype: 1}]},
						  {type: "Ketone", root: "C", bonds:[{el: "O", btype: 2}]},

						  {type: "Primary Amine", root: "N", bonds:[{el: "H", btype: 1},{el: "H", btype: 1},{el: "R", btype: 1}]},
						  {type: "Secondary Amine", root: "N", bonds:[{el: "H", btype: 1},{el: "R", btype: 1},{el: "R", btype: 1}]},
						  {type: "Tertiary Amine", root: "N", bonds:[{el: "R", btype: 1},{el: "R", btype: 1},{el: "R", btype: 1}]},

						  {type: "Nitro", root: "N", bonds:[{el: "O", btype: 2},{el: "O", btype: 1}]},
						  {type: "Alcohol", root: "O", bonds:[{el: "H", btype: 1}]},
						  {type: "Halo", root: "R", bonds:[{el: "X", btype: 1}]},
						  {type: "Nitrile", root: "C", bonds:[{el: "N", btype: 3}]},

						  {type: "Acetal", root: "C", bonds:[{el: "O", btype: 1}, {el: "O", btype: 1}]},
						  {type: "Ether", root: "O", bonds:[{el: "R", btype: 1},{el: "R", btype: 1}]},
						 ];
	exports.Atom = Atom;
	exports.Bond = Bond;
	exports.fGroup = fGroup;

}))
