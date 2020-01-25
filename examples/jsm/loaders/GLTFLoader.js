/**
 * @author Rich Tibbett / https://github.com/richtr
 * @author mrdoob / http://mrdoob.com/
 * @author Tony Parisi / http://www.tonyparisi.com/
 * @author Takahiro / https://github.com/takahirox
 * @author Don McCurdy / https://www.donmccurdy.com
 */

import {
	AnimationClip,
	Bone,
	Box3,
	BufferAttribute,
	BufferGeometry,
	ClampToEdgeWrapping,
	Color,
	DirectionalLight,
	DoubleSide,
	FileLoader,
	FrontSide,
	Group,
	InterleavedBuffer,
	InterleavedBufferAttribute,
	Interpolant,
	InterpolateDiscrete,
	InterpolateLinear,
	Line,
	LineBasicMaterial,
	LineLoop,
	LineSegments,
	LinearFilter,
	LinearMipmapLinearFilter,
	LinearMipmapNearestFilter,
	Loader,
	LoaderUtils,
	Material,
	MathUtils,
	Matrix4,
	Mesh,
	MeshBasicMaterial,
	MeshStandardMaterial,
	MirroredRepeatWrapping,
	NearestFilter,
	NearestMipmapLinearFilter,
	NearestMipmapNearestFilter,
	NumberKeyframeTrack,
	Object3D,
	OrthographicCamera,
	PerspectiveCamera,
	PointLight,
	Points,
	PointsMaterial,
	PropertyBinding,
	QuaternionKeyframeTrack,
	RGBAFormat,
	RGBFormat,
	RepeatWrapping,
	Scene,
	Skeleton,
	SkinnedMesh,
	Sphere,
	SpotLight,
	TextureLoader,
	TriangleFanDrawMode,
	TriangleStripDrawMode,
	Vector2,
	Vector3,
	VectorKeyframeTrack,
	VertexColors,
	sRGBEncoding
} from "../../../build/three.module.js";

var GLTFLoader = ( function () {

	function GLTFLoader( manager ) {

		Loader.call( this, manager );

		this.plugins = {
			extensions: {},
			generics: []
		};

		this.registerPlugin( new GLTFMaterialsPbrSpecularGlossinessExtension() );
		this.registerPlugin( new GLTFMaterialsUnlitExtension() );
		this.registerPlugin( new GLTFMeshQuantizationExtension() );
		this.registerPlugin( new GLTFLightsExtension() );
		this.registerPlugin( new GLTFTextureTransformExtension() );

	}

	GLTFLoader.prototype = Object.assign( Object.create( Loader.prototype ), {

		constructor: GLTFLoader,

		load: function ( url, onLoad, onProgress, onError ) {

			var scope = this;

			var resourcePath;

			if ( this.resourcePath !== '' ) {

				resourcePath = this.resourcePath;

			} else if ( this.path !== '' ) {

				resourcePath = this.path;

			} else {

				resourcePath = LoaderUtils.extractUrlBase( url );

			}

			// Tells the LoadingManager to track an extra item, which resolves after
			// the model is fully loaded. This means the count of items loaded will
			// be incorrect, but ensures manager.onLoad() does not fire early.
			scope.manager.itemStart( url );

			var _onError = function ( e ) {

				if ( onError ) {

					onError( e );

				} else {

					console.error( e );

				}

				scope.manager.itemError( url );
				scope.manager.itemEnd( url );

			};

			var loader = new FileLoader( scope.manager );

			loader.setPath( this.path );
			loader.setResponseType( 'arraybuffer' );

			if ( scope.crossOrigin === 'use-credentials' ) {

				loader.setWithCredentials( true );

			}

			loader.load( url, function ( data ) {

				try {

					scope.parse( data, resourcePath, function ( gltf ) {

						onLoad( gltf );

						scope.manager.itemEnd( url );

					}, _onError );

				} catch ( e ) {

					_onError( e );

				}

			}, onProgress, _onError );

		},

		setDRACOLoader: function ( dracoLoader ) {

			this.registerPlugin( new GLTFDracoMeshCompressionExtension( dracoLoader ) );
			return this;

		},

		registerPlugin: function ( plugin ) {

			var extensionName = plugin.extension;

			if ( extensionName ) {

				this.plugins.extensions[ extensionName ] = plugin;

			} else {

				if ( ! this.plugins.generics.includes( plugin ) ) {

					this.plugins.generics.push( plugin );

				}

			}

			return this;

		},

		setDDSLoader: function ( ddsLoader ) {

			this.registerPlugin( new GLTFTextureDDSExtension( ddsLoader ) );
			return this;

		},

		parse: function ( data, path, onLoad, onError ) {

			var content;
			var extensions = {};

			if ( typeof data === 'string' ) {

				content = data;

			} else {

				var magic = LoaderUtils.decodeText( new Uint8Array( data, 0, 4 ) );

				if ( magic === BINARY_EXTENSION_HEADER_MAGIC ) {

					try {

						// TODO: Binary is in core spec, not extension.
						// Maybe better to rename and to handle out of plugin system.
						this.registerPlugin( new GLTFBinaryExtension( data ) );

					} catch ( error ) {

						if ( onError ) onError( error );
						return;

					}

					content = this.plugins.extensions[ EXTENSIONS.KHR_BINARY_GLTF ].content;

				} else {

					content = LoaderUtils.decodeText( new Uint8Array( data ) );

				}

			}

			var json = JSON.parse( content );

			if ( json.asset === undefined || json.asset.version[ 0 ] < 2 ) {

				if ( onError ) onError( new Error( 'THREE.GLTFLoader: Unsupported asset. glTF versions >=2.0 are supported.' ) );
				return;

			}

			if ( json.extensionsUsed ) {

				for ( var i = 0; i < json.extensionsUsed.length; ++ i ) {

					var extensionName = json.extensionsUsed[ i ];
					var extensionsRequired = json.extensionsRequired || [];

					if ( extensionsRequired.indexOf( extensionName ) >= 0 && this.plugins.extensions[ extensionName ] === undefined ) {

						console.warn( 'THREE.GLTFLoader: Unknown extension "' + extensionName + '".' );

					}

				}

			}

			var parser = new GLTFParser( json, this.plugins, {

				path: path || this.resourcePath || '',
				crossOrigin: this.crossOrigin,
				manager: this.manager

			} );

			parser.parse( onLoad, onError );

		}

	} );

	/* GLTFREGISTRY */

	function GLTFRegistry() {

		var objects = {};

		return	{

			get: function ( key ) {

				return objects[ key ];

			},

			add: function ( key, object ) {

				objects[ key ] = object;

			},

			remove: function ( key ) {

				delete objects[ key ];

			},

			removeAll: function () {

				objects = {};

			}

		};

	}

	/*********************************/
	/********** EXTENSIONS ***********/
	/*********************************/

	var EXTENSIONS = {
		KHR_BINARY_GLTF: 'KHR_binary_glTF',
		KHR_DRACO_MESH_COMPRESSION: 'KHR_draco_mesh_compression',
		KHR_LIGHTS_PUNCTUAL: 'KHR_lights_punctual',
		KHR_MATERIALS_PBR_SPECULAR_GLOSSINESS: 'KHR_materials_pbrSpecularGlossiness',
		KHR_MATERIALS_UNLIT: 'KHR_materials_unlit',
		KHR_TEXTURE_TRANSFORM: 'KHR_texture_transform',
		KHR_MESH_QUANTIZATION: 'KHR_mesh_quantization',
		MSFT_TEXTURE_DDS: 'MSFT_texture_dds'
	};

	/**
	 * DDS Texture Extension
	 *
	 * Specification: https://github.com/KhronosGroup/glTF/tree/master/extensions/2.0/Vendor/MSFT_texture_dds
	 *
	 */
	function GLTFTextureDDSExtension( ddsLoader ) {

		if ( ! ddsLoader ) {

			throw new Error( 'THREE.GLTFLoader: Attempting to load .dds texture without importing DDSLoader' );

		}

		this.extension = EXTENSIONS.MSFT_TEXTURE_DDS;
		this.ddsLoader = ddsLoader;

	}

	GLTFTextureDDSExtension.prototype = {

		constructor: GLTFTextureDDSExtension,

		onTexture: function ( textureDef, parser ) {

			var extensions = textureDef.extensions;
			var json = parser.json;
			var imageDef = json.images[ extensions[ this.extension ].source ];

			return parser.loadTextureFile( imageDef, this.ddsLoader ).then( function ( texture ) {

				texture.flipY = false;

				if ( textureDef.name !== undefined ) texture.name = textureDef.name;

				var samplers = json.samplers || {};
				var sampler = samplers[ textureDef.sampler ] || {};

				texture.magFilter = WEBGL_FILTERS[ sampler.magFilter ] || LinearFilter;
				texture.minFilter = WEBGL_FILTERS[ sampler.minFilter ] || LinearMipmapLinearFilter;
				texture.wrapS = WEBGL_WRAPPINGS[ sampler.wrapS ] || RepeatWrapping;
				texture.wrapT = WEBGL_WRAPPINGS[ sampler.wrapT ] || RepeatWrapping;

				return texture;

			} );

		}

	};

	/**
	 * Punctual Lights Extension
	 *
	 * Specification: https://github.com/KhronosGroup/glTF/tree/master/extensions/2.0/Khronos/KHR_lights_punctual
	 */
	function GLTFLightsExtension() {

		this.extension = EXTENSIONS.KHR_LIGHTS_PUNCTUAL;

	}

	GLTFLightsExtension.prototype = {

		constructor: GLTFLightsExtension,

		onNode: function ( nodeDef, parser ) {

			var extensions = nodeDef.extensions || {};
			var json = parser.json;
			var lightExtension = json.extensions[ this.extension ] || {};
			var lightDefs = lightExtension.lights || [];
			var lightIndex = extensions[ this.extension ].light;
			var lightDef = lightDefs[ lightIndex ];

			var lightNode;

			var color = new Color( 0xffffff );

			if ( lightDef.color !== undefined ) color.fromArray( lightDef.color );

			var range = lightDef.range !== undefined ? lightDef.range : 0;

			switch ( lightDef.type ) {

				case 'directional':
					lightNode = new DirectionalLight( color );
					lightNode.target.position.set( 0, 0, - 1 );
					lightNode.add( lightNode.target );
					break;

				case 'point':
					lightNode = new PointLight( color );
					lightNode.distance = range;
					break;

				case 'spot':
					lightNode = new SpotLight( color );
					lightNode.distance = range;
					// Handle spotlight properties.
					lightDef.spot = lightDef.spot || {};
					lightDef.spot.innerConeAngle = lightDef.spot.innerConeAngle !== undefined ? lightDef.spot.innerConeAngle : 0;
					lightDef.spot.outerConeAngle = lightDef.spot.outerConeAngle !== undefined ? lightDef.spot.outerConeAngle : Math.PI / 4.0;
					lightNode.angle = lightDef.spot.outerConeAngle;
					lightNode.penumbra = 1.0 - lightDef.spot.innerConeAngle / lightDef.spot.outerConeAngle;
					lightNode.target.position.set( 0, 0, - 1 );
					lightNode.add( lightNode.target );
					break;

				default:
					throw new Error( 'THREE.GLTFLoader: Unexpected light type, "' + lightDef.type + '".' );

			}

			// Some lights (e.g. spot) default to a position other than the origin. Reset the position
			// here, because node-level parsing will only override position if explicitly specified.
			lightNode.position.set( 0, 0, 0 );

			lightNode.decay = 2;

			if ( lightDef.intensity !== undefined ) lightNode.intensity = lightDef.intensity;

			lightNode.name = lightDef.name || ( 'light_' + lightIndex );

			return Promise.resolve( lightNode );

		}

	};

	/**
	 * Unlit Materials Extension
	 *
	 * Specification: https://github.com/KhronosGroup/glTF/tree/master/extensions/2.0/Khronos/KHR_materials_unlit
	 */
	function GLTFMaterialsUnlitExtension() {

		this.extension = EXTENSIONS.KHR_MATERIALS_UNLIT;

	}

	GLTFMaterialsUnlitExtension.prototype = {

		constructor: GLTFMaterialsUnlitExtension,

		onMaterial: function ( materialDef, parser ) {

			var materialParams = {};
			var pending = [];

			materialParams.color = new Color( 1.0, 1.0, 1.0 );
			materialParams.opacity = 1.0;

			var metallicRoughness = materialDef.pbrMetallicRoughness;

			if ( metallicRoughness ) {

				if ( Array.isArray( metallicRoughness.baseColorFactor ) ) {

					var array = metallicRoughness.baseColorFactor;

					materialParams.color.fromArray( array );
					materialParams.opacity = array[ 3 ];

				}

				if ( metallicRoughness.baseColorTexture !== undefined ) {

					pending.push( parser.assignTexture( materialParams, 'map', metallicRoughness.baseColorTexture ) );

				}

			}

			return Promise.all( pending ).then( function () {

				return new MeshBasicMaterial( materialParams );

			} );

		}

	};

	/* BINARY EXTENSION */
	var BINARY_EXTENSION_HEADER_MAGIC = 'glTF';
	var BINARY_EXTENSION_HEADER_LENGTH = 12;
	var BINARY_EXTENSION_CHUNK_TYPES = { JSON: 0x4E4F534A, BIN: 0x004E4942 };

	function GLTFBinaryExtension( data ) {

		this.extension = EXTENSIONS.KHR_BINARY_GLTF;
		this.content = null;
		this.body = null;

		var headerView = new DataView( data, 0, BINARY_EXTENSION_HEADER_LENGTH );

		this.header = {
			magic: LoaderUtils.decodeText( new Uint8Array( data.slice( 0, 4 ) ) ),
			version: headerView.getUint32( 4, true ),
			length: headerView.getUint32( 8, true )
		};

		if ( this.header.magic !== BINARY_EXTENSION_HEADER_MAGIC ) {

			throw new Error( 'THREE.GLTFLoader: Unsupported glTF-Binary header.' );

		} else if ( this.header.version < 2.0 ) {

			throw new Error( 'THREE.GLTFLoader: Legacy binary file detected.' );

		}

		var chunkView = new DataView( data, BINARY_EXTENSION_HEADER_LENGTH );
		var chunkIndex = 0;

		while ( chunkIndex < chunkView.byteLength ) {

			var chunkLength = chunkView.getUint32( chunkIndex, true );
			chunkIndex += 4;

			var chunkType = chunkView.getUint32( chunkIndex, true );
			chunkIndex += 4;

			if ( chunkType === BINARY_EXTENSION_CHUNK_TYPES.JSON ) {

				var contentArray = new Uint8Array( data, BINARY_EXTENSION_HEADER_LENGTH + chunkIndex, chunkLength );
				this.content = LoaderUtils.decodeText( contentArray );

			} else if ( chunkType === BINARY_EXTENSION_CHUNK_TYPES.BIN ) {

				var byteOffset = BINARY_EXTENSION_HEADER_LENGTH + chunkIndex;
				this.body = data.slice( byteOffset, byteOffset + chunkLength );

			}

			// Clients must ignore chunks with unknown types.

			chunkIndex += chunkLength;

		}

		if ( this.content === null ) {

			throw new Error( 'THREE.GLTFLoader: JSON content not found.' );

		}

	}

	/**
	 * DRACO Mesh Compression Extension
	 *
	 * Specification: https://github.com/KhronosGroup/glTF/tree/master/extensions/2.0/Khronos/KHR_draco_mesh_compression
	 */
	function GLTFDracoMeshCompressionExtension( dracoLoader ) {

		if ( ! dracoLoader ) {

			throw new Error( 'THREE.GLTFLoader: No DRACOLoader instance provided.' );

		}

		this.extension = EXTENSIONS.KHR_DRACO_MESH_COMPRESSION;
		this.dracoLoader = dracoLoader;
		this.dracoLoader.preload();

	}

	GLTFDracoMeshCompressionExtension.prototype = {

		constructor: GLTFDracoMeshCompressionExtension,

		onGeometry: function ( primitive, parser ) {

			var json = parser.json;
			var dracoLoader = this.dracoLoader;
			var bufferViewIndex = primitive.extensions[ this.extension ].bufferView;
			var gltfAttributeMap = primitive.extensions[ this.extension ].attributes;
			var threeAttributeMap = {};
			var attributeNormalizedMap = {};
			var attributeTypeMap = {};

			for ( var attributeName in gltfAttributeMap ) {

				var threeAttributeName = ATTRIBUTES[ attributeName ] || attributeName.toLowerCase();

				threeAttributeMap[ threeAttributeName ] = gltfAttributeMap[ attributeName ];

			}

			for ( attributeName in primitive.attributes ) {

				var threeAttributeName = ATTRIBUTES[ attributeName ] || attributeName.toLowerCase();

				if ( gltfAttributeMap[ attributeName ] !== undefined ) {

					var accessorDef = json.accessors[ primitive.attributes[ attributeName ] ];
					var componentType = WEBGL_COMPONENT_TYPES[ accessorDef.componentType ];

					attributeTypeMap[ threeAttributeName ] = componentType;
					attributeNormalizedMap[ threeAttributeName ] = accessorDef.normalized === true;

				}

			}

			return parser.getDependency( 'bufferView', bufferViewIndex ).then( function ( bufferView ) {

				return new Promise( function ( resolve ) {

					dracoLoader.decodeDracoFile( bufferView, function ( geometry ) {

						for ( var attributeName in geometry.attributes ) {

							var attribute = geometry.attributes[ attributeName ];
							var normalized = attributeNormalizedMap[ attributeName ];

							if ( normalized !== undefined ) attribute.normalized = normalized;

						}

						resolve( geometry );

					}, threeAttributeMap, attributeTypeMap );

				} ).then( function ( geometry ) {

					return addPrimitiveAttributes( geometry, primitive, parser );

				} );

			} );

		},

	};

	/**
	 * Texture Transform Extension
	 *
	 * Specification: https://github.com/KhronosGroup/glTF/tree/master/extensions/2.0/Khronos/KHR_texture_transform
	 */
	function GLTFTextureTransformExtension() {

		this.extension = EXTENSIONS.KHR_TEXTURE_TRANSFORM;

	}

	GLTFTextureTransformExtension.prototype = {

		constructor: GLTFTextureTransformExtension,

		onAfterMap: function ( texture, textureDef, parser ) {

			var extensions = textureDef.extensions;
			var transform = extensions[ this.name ];

			texture = texture.clone();

			if ( transform.offset !== undefined ) {

				texture.offset.fromArray( transform.offset );

			}

			if ( transform.rotation !== undefined ) {

				texture.rotation = transform.rotation;

			}

			if ( transform.scale !== undefined ) {

				texture.repeat.fromArray( transform.scale );

			}

			if ( transform.texCoord !== undefined ) {

				console.warn( 'THREE.GLTFLoader: Custom UV sets in "' + this.name + '" extension not yet supported.' );

			}

			texture.needsUpdate = true;

			return texture;

		}

	};

	function GLTFMeshStandardSGMaterial( params ) {

		MeshStandardMaterial.call( this );

		var specularMapParsFragmentChunk = [
			'#ifdef USE_SPECULARMAP',
			'	uniform sampler2D specularMap;',
			'#endif'
		].join( '\n' );

		var glossinessMapParsFragmentChunk = [
			'#ifdef USE_GLOSSINESSMAP',
			'	uniform sampler2D glossinessMap;',
			'#endif'
		].join( '\n' );

		var specularMapFragmentChunk = [
			'vec3 specularFactor = specular;',
			'#ifdef USE_SPECULARMAP',
			'	vec4 texelSpecular = texture2D( specularMap, vUv );',
			'	texelSpecular = sRGBToLinear( texelSpecular );',
			'	// reads channel RGB, compatible with a glTF Specular-Glossiness (RGBA) texture',
			'	specularFactor *= texelSpecular.rgb;',
			'#endif'
		].join( '\n' );

		var glossinessMapFragmentChunk = [
			'float glossinessFactor = glossiness;',
			'#ifdef USE_GLOSSINESSMAP',
			'	vec4 texelGlossiness = texture2D( glossinessMap, vUv );',
			'	// reads channel A, compatible with a glTF Specular-Glossiness (RGBA) texture',
			'	glossinessFactor *= texelGlossiness.a;',
			'#endif'
		].join( '\n' );

		var lightPhysicalFragmentChunk = [
			'PhysicalMaterial material;',
			'material.diffuseColor = diffuseColor.rgb;',
			'material.specularRoughness = clamp( 1.0 - glossinessFactor, 0.04, 1.0 );',
			'material.specularColor = specularFactor.rgb;',
		].join( '\n' );

		var uniforms = {
			specular: { value: new Color().setHex( 0xffffff ) },
			glossiness: { value: 1 },
			specularMap: { value: null },
			glossinessMap: { value: null }
		};

		this._extraUniforms = uniforms;

		this.onBeforeCompile = function ( shader ) {

			for ( var uniformName in uniforms ) {

				shader.uniforms[ uniformName ] = uniforms[ uniformName ];

			}

			shader.fragmentShader = shader.fragmentShader.replace( 'uniform float roughness;', 'uniform vec3 specular;' );
			shader.fragmentShader = shader.fragmentShader.replace( 'uniform float metalness;', 'uniform float glossiness;' );
			shader.fragmentShader = shader.fragmentShader.replace( '#include <roughnessmap_pars_fragment>', specularMapParsFragmentChunk );
			shader.fragmentShader = shader.fragmentShader.replace( '#include <metalnessmap_pars_fragment>', glossinessMapParsFragmentChunk );
			shader.fragmentShader = shader.fragmentShader.replace( '#include <roughnessmap_fragment>', specularMapFragmentChunk );
			shader.fragmentShader = shader.fragmentShader.replace( '#include <metalnessmap_fragment>', glossinessMapFragmentChunk );
			shader.fragmentShader = shader.fragmentShader.replace( '#include <lights_physical_fragment>', lightPhysicalFragmentChunk );

		};

		delete this.metalness;
		delete this.roughness;
		delete this.metalnessMap;
		delete this.roughnessMap;

		this.setValues( params );

	}

	GLTFMeshStandardSGMaterial.prototype = Object.assign( Object.create( MeshStandardMaterial.prototype ), {

		constructor: GLTFMeshStandardSGMaterial,

		isGLTFMeshStandardSGMaterial: true,

		copy: function ( source ) {

			MeshStandardMaterial.prototype.copy.call( this, source );

			this.specularMap = source.specularMap;
			this.specular.copy( source.specular );
			this.glossinessMap = source.glossinessMap;
			this.glossiness = source.glossiness;

			delete this.metalness;
			delete this.roughness;
			delete this.metalnessMap;
			delete this.roughnessMap;

			return this;

		}

	} );

	Object.defineProperties( GLTFMeshStandardSGMaterial.prototype, {

		specular: {

			get: function () {

				return this._extraUniforms.specular.value;

			},

			set: function ( v ) {

				this._extraUniforms.specular.value = v;

			}

		},

		glossiness: {

			get: function () {

				return this._extraUniforms.glossiness.value;

			},

			set: function ( v ) {

				this._extraUniforms.glossiness.value = v;

			}

		},

		glossinessMap: {

			get: function () {

				this._extraUniforms.glossinessMap.value;

			},

			set: function ( v ) {

				this._extraUniforms.glossinessMap.value = v;

				if ( v ) {

					this.defines.USE_GLOSSINESSMAP = '';
					// set USE_ROUGHNESSMAP to enable vUv
					this.defines.USE_ROUGHNESSMAP = '';

				} else {

					delete this.defines.USE_ROUGHNESSMAP;
					delete this.defines.USE_GLOSSINESSMAP;

				}

			}

		}

	} );

	/**
	 * Specular-Glossiness Extension
	 *
	 * Specification: https://github.com/KhronosGroup/glTF/tree/master/extensions/2.0/Khronos/KHR_materials_pbrSpecularGlossiness
	 */
	function GLTFMaterialsPbrSpecularGlossinessExtension() {

		this.extension = EXTENSIONS.KHR_MATERIALS_PBR_SPECULAR_GLOSSINESS;

	}

	GLTFMaterialsPbrSpecularGlossinessExtension.prototype = {

		constructor: GLTFMaterialsPbrSpecularGlossinessExtension,

		onAfterMaterial: function( material, materialDef, parser ) {

			// Create Specular-Glossiness material from MeshStandardMaterial instance.
			// Rewrite to onMaterial if you want to avoid MeshStandardMaterial creation cost.

			var extensions = materialDef.extensions;
			var pbrSpecularGlossiness = extensions[ this.extension ];

			var newMaterial = new GLTFMeshStandardSGMaterial();

			MeshStandardMaterial.prototype.copy.call( newMaterial, material );

			delete newMaterial.metalness;
			delete newMaterial.roughness;
			delete newMaterial.metalnessMap;
			delete newMaterial.roughnessMap;

			var pending = [];

			newMaterial.color.setRGB( 1.0, 1.0, 1.0 );
			newMaterial.opacity = 1.0;

			if ( Array.isArray( pbrSpecularGlossiness.diffuseFactor ) ) {

				var array = pbrSpecularGlossiness.diffuseFactor;

				newMaterial.color.fromArray( array );
				newMaterial.opacity = array[ 3 ];

			}

			newMaterial.map = null;

			if ( pbrSpecularGlossiness.diffuseTexture !== undefined ) {

				pending.push( parser.assignTexture( newMaterial, 'map', pbrSpecularGlossiness.diffuseTexture ) );

			}

			newMaterial.glossiness = pbrSpecularGlossiness.glossinessFactor !== undefined ? pbrSpecularGlossiness.glossinessFactor : 1.0;
			newMaterial.specular.setRGB( 1.0, 1.0, 1.0 );

			if ( Array.isArray( pbrSpecularGlossiness.specularFactor ) ) {

				newMaterial.specular.fromArray( pbrSpecularGlossiness.specularFactor );

			}

			if ( pbrSpecularGlossiness.specularGlossinessTexture !== undefined ) {

				var specGlossMapDef = pbrSpecularGlossiness.specularGlossinessTexture;
				pending.push( parser.assignTexture( newMaterial, 'glossinessMap', specGlossMapDef ) );
				pending.push( parser.assignTexture( newMaterial, 'specularMap', specGlossMapDef ) );

			}

			return Promise.all( pending ).then( function () {

				if ( newMaterial.map ) newMaterial.map.encoding = sRGBEncoding;
				if ( newMaterial.specularMap ) newMaterial.specularMap.encoding = sRGBEncoding;

				return newMaterial;

			} );

		}

	};

	/**
	 * Mesh Quantization Extension
	 *
	 * Specification: https://github.com/KhronosGroup/glTF/tree/master/extensions/2.0/Khronos/KHR_mesh_quantization
	 */
	function GLTFMeshQuantizationExtension() {

		this.extension = EXTENSIONS.KHR_MESH_QUANTIZATION;

	}

	/*********************************/
	/********** INTERPOLATION ********/
	/*********************************/

	// Spline Interpolation
	// Specification: https://github.com/KhronosGroup/glTF/blob/master/specification/2.0/README.md#appendix-c-spline-interpolation
	function GLTFCubicSplineInterpolant( parameterPositions, sampleValues, sampleSize, resultBuffer ) {

		Interpolant.call( this, parameterPositions, sampleValues, sampleSize, resultBuffer );

	}

	GLTFCubicSplineInterpolant.prototype = Object.create( Interpolant.prototype );
	GLTFCubicSplineInterpolant.prototype.constructor = GLTFCubicSplineInterpolant;

	GLTFCubicSplineInterpolant.prototype.copySampleValue_ = function ( index ) {

		// Copies a sample value to the result buffer. See description of glTF
		// CUBICSPLINE values layout in interpolate_() function below.

		var result = this.resultBuffer,
			values = this.sampleValues,
			valueSize = this.valueSize,
			offset = index * valueSize * 3 + valueSize;

		for ( var i = 0; i !== valueSize; i ++ ) {

			result[ i ] = values[ offset + i ];

		}

		return result;

	};

	GLTFCubicSplineInterpolant.prototype.beforeStart_ = GLTFCubicSplineInterpolant.prototype.copySampleValue_;

	GLTFCubicSplineInterpolant.prototype.afterEnd_ = GLTFCubicSplineInterpolant.prototype.copySampleValue_;

	GLTFCubicSplineInterpolant.prototype.interpolate_ = function ( i1, t0, t, t1 ) {

		var result = this.resultBuffer;
		var values = this.sampleValues;
		var stride = this.valueSize;

		var stride2 = stride * 2;
		var stride3 = stride * 3;

		var td = t1 - t0;

		var p = ( t - t0 ) / td;
		var pp = p * p;
		var ppp = pp * p;

		var offset1 = i1 * stride3;
		var offset0 = offset1 - stride3;

		var s2 = - 2 * ppp + 3 * pp;
		var s3 = ppp - pp;
		var s0 = 1 - s2;
		var s1 = s3 - pp + p;

		// Layout of keyframe output values for CUBICSPLINE animations:
		//   [ inTangent_1, splineVertex_1, outTangent_1, inTangent_2, splineVertex_2, ... ]
		for ( var i = 0; i !== stride; i ++ ) {

			var p0 = values[ offset0 + i + stride ]; // splineVertex_k
			var m0 = values[ offset0 + i + stride2 ] * td; // outTangent_k * (t_k+1 - t_k)
			var p1 = values[ offset1 + i + stride ]; // splineVertex_k+1
			var m1 = values[ offset1 + i ] * td; // inTangent_k+1 * (t_k+1 - t_k)

			result[ i ] = s0 * p0 + s1 * m0 + s2 * p1 + s3 * m1;

		}

		return result;

	};

	/*********************************/
	/********** INTERNALS ************/
	/*********************************/

	/* CONSTANTS */

	var WEBGL_CONSTANTS = {
		FLOAT: 5126,
		//FLOAT_MAT2: 35674,
		FLOAT_MAT3: 35675,
		FLOAT_MAT4: 35676,
		FLOAT_VEC2: 35664,
		FLOAT_VEC3: 35665,
		FLOAT_VEC4: 35666,
		LINEAR: 9729,
		REPEAT: 10497,
		SAMPLER_2D: 35678,
		POINTS: 0,
		LINES: 1,
		LINE_LOOP: 2,
		LINE_STRIP: 3,
		TRIANGLES: 4,
		TRIANGLE_STRIP: 5,
		TRIANGLE_FAN: 6,
		UNSIGNED_BYTE: 5121,
		UNSIGNED_SHORT: 5123
	};

	var WEBGL_COMPONENT_TYPES = {
		5120: Int8Array,
		5121: Uint8Array,
		5122: Int16Array,
		5123: Uint16Array,
		5125: Uint32Array,
		5126: Float32Array
	};

	var WEBGL_FILTERS = {
		9728: NearestFilter,
		9729: LinearFilter,
		9984: NearestMipmapNearestFilter,
		9985: LinearMipmapNearestFilter,
		9986: NearestMipmapLinearFilter,
		9987: LinearMipmapLinearFilter
	};

	var WEBGL_WRAPPINGS = {
		33071: ClampToEdgeWrapping,
		33648: MirroredRepeatWrapping,
		10497: RepeatWrapping
	};

	var WEBGL_TYPE_SIZES = {
		'SCALAR': 1,
		'VEC2': 2,
		'VEC3': 3,
		'VEC4': 4,
		'MAT2': 4,
		'MAT3': 9,
		'MAT4': 16
	};

	var ATTRIBUTES = {
		POSITION: 'position',
		NORMAL: 'normal',
		TANGENT: 'tangent',
		TEXCOORD_0: 'uv',
		TEXCOORD_1: 'uv2',
		COLOR_0: 'color',
		WEIGHTS_0: 'skinWeight',
		JOINTS_0: 'skinIndex',
	};

	var PATH_PROPERTIES = {
		scale: 'scale',
		translation: 'position',
		rotation: 'quaternion',
		weights: 'morphTargetInfluences'
	};

	var INTERPOLATION = {
		CUBICSPLINE: undefined, // We use a custom interpolant (GLTFCubicSplineInterpolation) for CUBICSPLINE tracks. Each
		                        // keyframe track will be initialized with a default interpolation type, then modified.
		LINEAR: InterpolateLinear,
		STEP: InterpolateDiscrete
	};

	var ALPHA_MODES = {
		OPAQUE: 'OPAQUE',
		MASK: 'MASK',
		BLEND: 'BLEND'
	};

	var MIME_TYPE_FORMATS = {
		'image/png': RGBAFormat,
		'image/jpeg': RGBFormat
	};

	/* UTILITY FUNCTIONS */

	function resolveURL( url, path ) {

		// Invalid URL
		if ( typeof url !== 'string' || url === '' ) return '';

		// Host Relative URL
		if ( /^https?:\/\//i.test( path ) && /^\//.test( url ) ) {

			path = path.replace( /(^https?:\/\/[^\/]+).*/i, '$1' );

		}

		// Absolute URL http://,https://,//
		if ( /^(https?:)?\/\//i.test( url ) ) return url;

		// Data URI
		if ( /^data:.*,.*$/i.test( url ) ) return url;

		// Blob URL
		if ( /^blob:.*$/i.test( url ) ) return url;

		// Relative URL
		return path + url;

	}

	/**
	 * Specification: https://github.com/KhronosGroup/glTF/blob/master/specification/2.0/README.md#default-material
	 */
	function createDefaultMaterial( cache ) {

		if ( cache[ 'DefaultMaterial' ] === undefined ) {

			cache[ 'DefaultMaterial' ] = new MeshStandardMaterial( {
				color: 0xFFFFFF,
				emissive: 0x000000,
				metalness: 1,
				roughness: 1,
				transparent: false,
				depthTest: true,
				side: FrontSide
			} );

		}

		return cache[ 'DefaultMaterial' ];

	}

	function addUnknownExtensionsToUserData( knownExtensions, object, objectDef ) {

		// Add unknown glTF extensions to an object's userData.

		for ( var name in objectDef.extensions ) {

			if ( knownExtensions[ name ] === undefined ) {

				object.userData.gltfExtensions = object.userData.gltfExtensions || {};
				object.userData.gltfExtensions[ name ] = objectDef.extensions[ name ];

			}

		}

	}

	/**
	 * @param {Object3D|Material|BufferGeometry} object
	 * @param {GLTF.definition} gltfDef
	 */
	function assignExtrasToUserData( object, gltfDef ) {

		if ( gltfDef.extras !== undefined ) {

			if ( typeof gltfDef.extras === 'object' ) {

				Object.assign( object.userData, gltfDef.extras );

			} else {

				console.warn( 'THREE.GLTFLoader: Ignoring primitive type .extras, ' + gltfDef.extras );

			}

		}

	}

	/**
	 * Specification: https://github.com/KhronosGroup/glTF/blob/master/specification/2.0/README.md#morph-targets
	 *
	 * @param {BufferGeometry} geometry
	 * @param {Array<GLTF.Target>} targets
	 * @param {GLTFParser} parser
	 * @return {Promise<BufferGeometry>}
	 */
	function addMorphTargets( geometry, targets, parser ) {

		var hasMorphPosition = false;
		var hasMorphNormal = false;

		for ( var i = 0, il = targets.length; i < il; i ++ ) {

			var target = targets[ i ];

			if ( target.POSITION !== undefined ) hasMorphPosition = true;
			if ( target.NORMAL !== undefined ) hasMorphNormal = true;

			if ( hasMorphPosition && hasMorphNormal ) break;

		}

		if ( ! hasMorphPosition && ! hasMorphNormal ) return Promise.resolve( geometry );

		var pendingPositionAccessors = [];
		var pendingNormalAccessors = [];

		for ( var i = 0, il = targets.length; i < il; i ++ ) {

			var target = targets[ i ];

			if ( hasMorphPosition ) {

				var pendingAccessor = target.POSITION !== undefined
					? parser.getDependency( 'accessor', target.POSITION )
					: geometry.attributes.position;

				pendingPositionAccessors.push( pendingAccessor );

			}

			if ( hasMorphNormal ) {

				var pendingAccessor = target.NORMAL !== undefined
					? parser.getDependency( 'accessor', target.NORMAL )
					: geometry.attributes.normal;

				pendingNormalAccessors.push( pendingAccessor );

			}

		}

		return Promise.all( [
			Promise.all( pendingPositionAccessors ),
			Promise.all( pendingNormalAccessors )
		] ).then( function ( accessors ) {

			var morphPositions = accessors[ 0 ];
			var morphNormals = accessors[ 1 ];

			if ( hasMorphPosition ) geometry.morphAttributes.position = morphPositions;
			if ( hasMorphNormal ) geometry.morphAttributes.normal = morphNormals;
			geometry.morphTargetsRelative = true;

			return geometry;

		} );

	}

	/**
	 * @param {Mesh} mesh
	 * @param {GLTF.Mesh} meshDef
	 */
	function updateMorphTargets( mesh, meshDef ) {

		mesh.updateMorphTargets();

		if ( meshDef.weights !== undefined ) {

			for ( var i = 0, il = meshDef.weights.length; i < il; i ++ ) {

				mesh.morphTargetInfluences[ i ] = meshDef.weights[ i ];

			}

		}

		// .extras has user-defined data, so check that .extras.targetNames is an array.
		if ( meshDef.extras && Array.isArray( meshDef.extras.targetNames ) ) {

			var targetNames = meshDef.extras.targetNames;

			if ( mesh.morphTargetInfluences.length === targetNames.length ) {

				mesh.morphTargetDictionary = {};

				for ( var i = 0, il = targetNames.length; i < il; i ++ ) {

					mesh.morphTargetDictionary[ targetNames[ i ] ] = i;

				}

			} else {

				console.warn( 'THREE.GLTFLoader: Invalid extras.targetNames length. Ignoring names.' );

			}

		}

	}

	function createPrimitiveKey( primitiveDef ) {

		var dracoExtension = primitiveDef.extensions && primitiveDef.extensions[ EXTENSIONS.KHR_DRACO_MESH_COMPRESSION ];
		var geometryKey;

		if ( dracoExtension ) {

			// TODO: Can we move this under GLTFDracoMeshCompressionExtension?
			geometryKey = 'draco:' + dracoExtension.bufferView
				+ ':' + dracoExtension.indices
				+ ':' + createAttributesKey( dracoExtension.attributes );

		} else {

			geometryKey = primitiveDef.indices + ':' + createAttributesKey( primitiveDef.attributes ) + ':' + primitiveDef.mode;

		}

		return geometryKey;

	}

	function createAttributesKey( attributes ) {

		var attributesKey = '';

		var keys = Object.keys( attributes ).sort();

		for ( var i = 0, il = keys.length; i < il; i ++ ) {

			attributesKey += keys[ i ] + ':' + attributes[ keys[ i ] ] + ';';

		}

		return attributesKey;

	}

	/* GLTF PARSER */

	function GLTFParser( json, plugins, options ) {

		this.json = json || {};
		this.plugins = plugins || { extensions: {}, generics: [] };
		this.options = options || {};

		// loader object cache
		this.cache = new GLTFRegistry();

		// BufferGeometry caching
		this.primitiveCache = {};

		this.textureLoader = new TextureLoader( this.options.manager );
		this.textureLoader.setCrossOrigin( this.options.crossOrigin );

		this.fileLoader = new FileLoader( this.options.manager );
		this.fileLoader.setResponseType( 'arraybuffer' );

		if ( this.options.crossOrigin === 'use-credentials' ) {

			this.fileLoader.setWithCredentials( true );

		}

	}

	GLTFParser.prototype.parse = function ( onLoad, onError ) {

		var parser = this;
		var json = this.json;
		var extensions = this.plugins.extensions;

		// Clear the loader cache
		this.cache.removeAll();

		// Mark the special nodes/meshes in json for efficient parse
		this.markDefs();

		this._onBefore( 'GLTF', json ).then( function ( json ) {

			parser.json = json;

			return Promise.all( [

				parser.getDependencies( 'scene' ),
				parser.getDependencies( 'animation' ),
				parser.getDependencies( 'camera' ),

			] );

		} ).then( function ( dependencies ) {

			var result = {
				scene: dependencies[ 0 ][ json.scene || 0 ],
				scenes: dependencies[ 0 ],
				animations: dependencies[ 1 ],
				cameras: dependencies[ 2 ],
				asset: json.asset,
				parser: parser,
				userData: {}
			};

			addUnknownExtensionsToUserData( extensions, result, json );

			assignExtrasToUserData( result, json );

			parser._onAfter( 'GLTF', result, json ).then( onLoad );

		} ).catch( onError );

	};

	/**
	 * Marks the special nodes/meshes in json for efficient parse.
	 */
	GLTFParser.prototype.markDefs = function () {

		var nodeDefs = this.json.nodes || [];
		var skinDefs = this.json.skins || [];
		var meshDefs = this.json.meshes || [];

		var meshReferences = {};
		var meshUses = {};

		// Nothing in the node definition indicates whether it is a Bone or an
		// Object3D. Use the skins' joint references to mark bones.
		for ( var skinIndex = 0, skinLength = skinDefs.length; skinIndex < skinLength; skinIndex ++ ) {

			var joints = skinDefs[ skinIndex ].joints;

			for ( var i = 0, il = joints.length; i < il; i ++ ) {

				nodeDefs[ joints[ i ] ].isBone = true;

			}

		}

		// Meshes can (and should) be reused by multiple nodes in a glTF asset. To
		// avoid having more than one Mesh with the same name, count
		// references and rename instances below.
		//
		// Example: CesiumMilkTruck sample model reuses "Wheel" meshes.
		for ( var nodeIndex = 0, nodeLength = nodeDefs.length; nodeIndex < nodeLength; nodeIndex ++ ) {

			var nodeDef = nodeDefs[ nodeIndex ];

			if ( nodeDef.mesh !== undefined ) {

				if ( meshReferences[ nodeDef.mesh ] === undefined ) {

					meshReferences[ nodeDef.mesh ] = meshUses[ nodeDef.mesh ] = 0;

				}

				meshReferences[ nodeDef.mesh ] ++;

				// Nothing in the mesh definition indicates whether it is
				// a SkinnedMesh or Mesh. Use the node's mesh reference
				// to mark SkinnedMesh if node has skin.
				if ( nodeDef.skin !== undefined ) {

					meshDefs[ nodeDef.mesh ].isSkinnedMesh = true;

				}

			}

		}

		this.json.meshReferences = meshReferences;
		this.json.meshUses = meshUses;

	};

	/**
	 * @param {string} key
	 * @param {GLTF.definition} def
	 * @return {Promise<GLTF.definition>}}
	 */
	GLTFParser.prototype._onBefore = function ( key, def ) {

		var parser = this;
		var functionName = 'onBefore' + key;
		var plugins = this.plugins || {};
		var extensionPlugins = plugins.extensions || {};
		var genericPlugins = plugins.generics || {};
		var extensions = def.extensions || {};
		var pending = Promise.resolve( def );

		for ( var extensionName in extensions ) {

			var plugin = extensionPlugins[ extensionName ];

			if ( plugin === undefined || plugin[ functionName ] === undefined ) {

				continue;

			}

			pending = pending.then( function ( plugin, def ) {

				return plugin[ functionName ]( def, parser );

			}.bind( null, plugin ) );

		}

		for ( var plugin of genericPlugins ) {

			if ( plugin[ functionName ] === undefined ) {

				continue;

			}

			pending = pending.then( function ( plugin, def ) {

				return plugin[ functionName ]( def, parser );

			}.bind( null, plugin ) );

		}

		return pending;

	};

	/**
	 * @param {string} key
	 * @param {Object} object
	 * @param {GLTF.definition} def
	 * @return {Promise<Object>}}
	 */
	GLTFParser.prototype._onAfter = function ( key, object, def ) {

		var parser = this;
		var functionName = 'onAfter' + key;
		var plugins = this.plugins || {};
		var extensionPlugins = plugins.extensions || {};
		var genericPlugins = plugins.generics || [];
		var extensions = def.extensions || {};
		var pending = Promise.resolve( object );

		for ( var extensionName in extensions ) {

			var plugin = extensionPlugins[ extensionName ];

			if ( plugin === undefined || plugin[ functionName ] === undefined ) {

				continue;

			}

			pending = pending.then( function ( plugin, object ) {

				return plugin[ functionName ]( object, def, parser );

			}.bind( null, plugin ) );

		}

		for ( var plugin of genericPlugins ) {

			if ( plugin[ functionName ] === undefined ) {

				continue;

			}

			pending = pending.then( function ( plugin, object ) {

				return plugin[ functionName ]( object, def, parser );

			}.bind( null, plugin ) );

		}

		return pending;

	};

	/**
	 * @param {string} key
	 * @param {GLTF.definition} def
	 * @return {Promise<Object>|null}}
	 */
	GLTFParser.prototype._on = function ( key, def ) {

		var parser = this;
		var functionName = 'on' + key;
		var plugins = this.plugins || {};
		var extensionPlugins = plugins.extensions || {};
		var genericPlugins = plugins.generics || [];
		var extensions = def.extensions || {};

		for ( var extensionName in extensions ) {

			var plugin = extensionPlugins[ extensionName ];

			if ( plugin === undefined || plugin[ functionName ] === undefined ) {

				continue;

			}

			var result = plugin[ functionName ]( def, parser );

			if ( result ) return result;

		}

		for ( var plugin of genericPlugins ) {

			if ( plugin[ functionName ] === undefined ) {

				continue;

			}

			var result = plugin[ functionName ]( def, parser );

			if ( result ) return result;

		}

		return Promise.resolve( null );

	};

	/**
	 * Requests the specified dependency asynchronously, with caching.
	 * @param {string} type
	 * @param {number} index
	 * @return {Promise<Object3D|Material|THREE.Texture|AnimationClip|ArrayBuffer|Object>}
	 */
	GLTFParser.prototype.getDependency = function ( type, index ) {

		var parser = this;
		var json = parser.json;
		var cacheKey = type + ':' + index;
		var dependency = this.cache.get( cacheKey );

		if ( ! dependency ) {

			switch ( type ) {

				case 'scene':
					dependency = this._onBefore( 'Scene', json.scenes[ index ] ).then( function ( def ) {

						return parser._on( 'Scene', def ).then( function ( customScene ) {

							return parser.loadScene( index, def, customScene );

						} ).then( function ( scene ) {

							return parser._onAfter( 'Scene', scene, def );

						} );

					} );
					break;

				case 'node':
					dependency = this._onBefore( 'Node', json.nodes[ index ] ).then( function ( def ) {

						return parser._on( 'Node', def ).then( function ( customNode ) {

							return parser.loadNode( index, def, customNode );

						} ).then( function ( node ) {

							return parser._onAfter( 'Node', node, def );

						} );

					} );
					break;

				case 'mesh':
					dependency = this._onBefore( 'Mesh', json.nodes[ index ] ).then( function ( def ) {

						return parser._on( 'Mesh', def ).then( function ( customMesh ) {

							return parser.loadMesh( index, def, customMesh );

						} ).then( function ( mesh ) {

							return parser._onAfter( 'Mesh', mesh, def );

						} );

					} );
					break;

				case 'accessor':
					dependency = this._onBefore( 'Accessor', json.accessors[ index ] ).then( function ( def ) {

						return parser._on( 'Accessor', def ).then( function ( accessor ) {

							return accessor || parser.loadAccessor( def );

						} ).then( function ( accessor ) {

							return parser._onAfter( 'Accessor', accessor, def );

						} );

					} );
					break;

				case 'bufferView':
					dependency = this._onBefore( 'BufferView', json.bufferViews[ index ] ).then( function ( def ) {

						return parser._on( 'BufferView', def ).then( function ( bufferView ) {

							return bufferView || parser.loadBufferView( def );

						} ).then( function ( bufferView ) {

							return parser._onAfter( 'BufferView', bufferView, def );

						} );

					} );
					break;

				case 'buffer':
					dependency = this._onBefore( 'Buffer', json.buffers[ index ] ).then( function ( def ) {

						return parser._on( 'Buffer', def ).then( function ( buffer ) {

							return buffer || parser.loadBuffer( index, def );

						} ).then( function ( buffer ) {

							return parser._onAfter( 'Buffer', buffer, def );

						} );

					} );
					break;

				case 'material':
					dependency = this._onBefore( 'Material', json.materials[ index ] ).then( function ( def ) {

						return parser._on( 'Material', def ).then( function ( customMaterial ) {

							return parser.loadMaterial( index, def, customMaterial );

						} ).then( function ( material ) {

							return parser._onAfter( 'Material', material, def );

						} );

					} );
					break;

				case 'texture':
					dependency = this._onBefore( 'Texture', json.textures[ index ] ).then( function ( def ) {

						return parser._on( 'Texture', def ).then( function ( texture ) {

							return texture || parser.loadTexture( def );

						} ).then( function ( texture ) {

							return parser._onAfter( 'Texture', texture, def );

						} );

					} );
					break;

				case 'skin':
					dependency = this._onBefore( 'Skin', json.skins[ index ] ).then( function ( def ) {

						return parser._on( 'Skin', def ).then( function ( skin ) {

							return skin || parser.loadSkin( def );

						} ).then( function ( skin ) {

							return parser._onAfter( 'Skin', skin, def );

						} );

					} );
					break;

				case 'animation':
					dependency = this._onBefore( 'Animation', json.animations[ index ] ).then( function ( def ) {

						return parser._on( 'Animation', def ).then( function ( animation ) {

							return animation || parser.loadAnimation( index, def );

						} ).then( function ( animation ) {

							return parser._onAfter( 'Animation', animation, def );

						} );

					} );
					break;

				case 'camera':
					dependency = this._onBefore( 'Camera', json.cameras[ index ] ).then( function ( def ) {

						return parser._on( 'Camera', def ).then( function ( camera ) {

							return camera || parser.loadCamera( def );

						} ).then( function ( camera ) {

							return parser._onAfter( 'Camera', camera, def );

						} );

					} );
					break;

				default:
					throw new Error( 'Unknown type: ' + type );

			}

			this.cache.add( cacheKey, dependency );

		}

		return dependency;

	};

	/**
	 * Requests all dependencies of the specified type asynchronously, with caching.
	 * @param {string} type
	 * @return {Promise<Array<Object>>}
	 */
	GLTFParser.prototype.getDependencies = function ( type ) {

		var dependencies = this.cache.get( type );

		if ( ! dependencies ) {

			var parser = this;
			var defs = this.json[ type + ( type === 'mesh' ? 'es' : 's' ) ] || [];

			dependencies = Promise.all( defs.map( function ( def, index ) {

				return parser.getDependency( type, index );

			} ) );

			this.cache.add( type, dependencies );

		}

		return dependencies;

	};

	/**
	 * Specification: https://github.com/KhronosGroup/glTF/blob/master/specification/2.0/README.md#buffers-and-buffer-views
	 * @param {number} bufferIndex
	 * @param {GLTF.definition} bufferDef
	 * @return {Promise<ArrayBuffer>}
	 */
	GLTFParser.prototype.loadBuffer = function ( bufferIndex, bufferDef ) {

		var loader = this.fileLoader;

		if ( bufferDef.type && bufferDef.type !== 'arraybuffer' ) {

			throw new Error( 'THREE.GLTFLoader: ' + bufferDef.type + ' buffer type is not supported.' );

		}

		// If present, GLB container is required to be the first buffer.
		if ( bufferDef.uri === undefined && bufferIndex === 0 ) {

			return Promise.resolve( this.plugins.extensions[ EXTENSIONS.KHR_BINARY_GLTF ].body );

		}

		var options = this.options;

		return new Promise( function ( resolve, reject ) {

			loader.load( resolveURL( bufferDef.uri, options.path ), resolve, undefined, function () {

				reject( new Error( 'THREE.GLTFLoader: Failed to load buffer "' + bufferDef.uri + '".' ) );

			} );

		} );

	};

	/**
	 * Specification: https://github.com/KhronosGroup/glTF/blob/master/specification/2.0/README.md#buffers-and-buffer-views
	 * @param {GLTF.definition} bufferViewDef
	 * @return {Promise<ArrayBuffer>}
	 */
	GLTFParser.prototype.loadBufferView = function ( bufferViewDef ) {

		return this.getDependency( 'buffer', bufferViewDef.buffer ).then( function ( buffer ) {

			var byteLength = bufferViewDef.byteLength || 0;
			var byteOffset = bufferViewDef.byteOffset || 0;
			return buffer.slice( byteOffset, byteOffset + byteLength );

		} );

	};

	/**
	 * Specification: https://github.com/KhronosGroup/glTF/blob/master/specification/2.0/README.md#accessors
	 * @param [GLTF.definition} accessorDef
	 * @return {Promise<BufferAttribute|InterleavedBufferAttribute>}
	 */
	GLTFParser.prototype.loadAccessor = function ( accessorDef ) {

		var parser = this;
		var json = this.json;

		if ( accessorDef.bufferView === undefined && accessorDef.sparse === undefined ) {

			// Ignore empty accessors, which may be used to declare runtime
			// information about attributes coming from another source (e.g. Draco
			// compression extension).
			return Promise.resolve( null );

		}

		var pendingBufferViews = [];

		if ( accessorDef.bufferView !== undefined ) {

			pendingBufferViews.push( this.getDependency( 'bufferView', accessorDef.bufferView ) );

		} else {

			pendingBufferViews.push( null );

		}

		if ( accessorDef.sparse !== undefined ) {

			pendingBufferViews.push( this.getDependency( 'bufferView', accessorDef.sparse.indices.bufferView ) );
			pendingBufferViews.push( this.getDependency( 'bufferView', accessorDef.sparse.values.bufferView ) );

		}

		return Promise.all( pendingBufferViews ).then( function ( bufferViews ) {

			var bufferView = bufferViews[ 0 ];

			var itemSize = WEBGL_TYPE_SIZES[ accessorDef.type ];
			var TypedArray = WEBGL_COMPONENT_TYPES[ accessorDef.componentType ];

			// For VEC3: itemSize is 3, elementBytes is 4, itemBytes is 12.
			var elementBytes = TypedArray.BYTES_PER_ELEMENT;
			var itemBytes = elementBytes * itemSize;
			var byteOffset = accessorDef.byteOffset || 0;
			var byteStride = accessorDef.bufferView !== undefined ? json.bufferViews[ accessorDef.bufferView ].byteStride : undefined;
			var normalized = accessorDef.normalized === true;
			var array, bufferAttribute;

			// The buffer is not interleaved if the stride is the item size in bytes.
			if ( byteStride && byteStride !== itemBytes ) {

				// Each "slice" of the buffer, as defined by 'count' elements of 'byteStride' bytes, gets its own InterleavedBuffer
				// This makes sure that IBA.count reflects accessor.count properly
				var ibSlice = Math.floor( byteOffset / byteStride );
				var ibCacheKey = 'InterleavedBuffer:' + accessorDef.bufferView + ':' + accessorDef.componentType + ':' + ibSlice + ':' + accessorDef.count;
				var ib = parser.cache.get( ibCacheKey );

				if ( ! ib ) {

					array = new TypedArray( bufferView, ibSlice * byteStride, accessorDef.count * byteStride / elementBytes );

					// Integer parameters to IB/IBA are in array elements, not bytes.
					ib = new InterleavedBuffer( array, byteStride / elementBytes );

					parser.cache.add( ibCacheKey, ib );

				}

				bufferAttribute = new InterleavedBufferAttribute( ib, itemSize, ( byteOffset % byteStride ) / elementBytes, normalized );

			} else {

				if ( bufferView === null ) {

					array = new TypedArray( accessorDef.count * itemSize );

				} else {

					array = new TypedArray( bufferView, byteOffset, accessorDef.count * itemSize );

				}

				bufferAttribute = new BufferAttribute( array, itemSize, normalized );

			}

			// https://github.com/KhronosGroup/glTF/blob/master/specification/2.0/README.md#sparse-accessors
			if ( accessorDef.sparse !== undefined ) {

				var itemSizeIndices = WEBGL_TYPE_SIZES.SCALAR;
				var TypedArrayIndices = WEBGL_COMPONENT_TYPES[ accessorDef.sparse.indices.componentType ];

				var byteOffsetIndices = accessorDef.sparse.indices.byteOffset || 0;
				var byteOffsetValues = accessorDef.sparse.values.byteOffset || 0;

				var sparseIndices = new TypedArrayIndices( bufferViews[ 1 ], byteOffsetIndices, accessorDef.sparse.count * itemSizeIndices );
				var sparseValues = new TypedArray( bufferViews[ 2 ], byteOffsetValues, accessorDef.sparse.count * itemSize );

				if ( bufferView !== null ) {

					// Avoid modifying the original ArrayBuffer, if the bufferView wasn't initialized with zeroes.
					bufferAttribute = new BufferAttribute( bufferAttribute.array.slice(), bufferAttribute.itemSize, bufferAttribute.normalized );

				}

				for ( var i = 0, il = sparseIndices.length; i < il; i ++ ) {

					var index = sparseIndices[ i ];

					bufferAttribute.setX( index, sparseValues[ i * itemSize ] );
					if ( itemSize >= 2 ) bufferAttribute.setY( index, sparseValues[ i * itemSize + 1 ] );
					if ( itemSize >= 3 ) bufferAttribute.setZ( index, sparseValues[ i * itemSize + 2 ] );
					if ( itemSize >= 4 ) bufferAttribute.setW( index, sparseValues[ i * itemSize + 3 ] );
					if ( itemSize >= 5 ) throw new Error( 'THREE.GLTFLoader: Unsupported itemSize in sparse BufferAttribute.' );

				}

			}

			return bufferAttribute;

		} );

	};

	/**
	 * @param {GLTF.Image} imageDef
	 * @param {Loader} loader (optional)
	 * @return {Promise<THREE.Texture>}
	 */
	GLTFParser.prototype.loadTextureFile = function ( imageDef, loader ) {

		var parser = this;
		var options = this.options;

		var sourceURI = imageDef.uri;
		var isObjectURL = false;

		if ( imageDef.bufferView !== undefined ) {

			// Load binary image data from bufferView, if provided.

			var URL = window.URL || window.webkitURL;

			sourceURI = parser.getDependency( 'bufferView', imageDef.bufferView ).then( function ( bufferView ) {

				isObjectURL = true;
				var blob = new Blob( [ bufferView ], { type: imageDef.mimeType } );
				return URL.createObjectURL( blob );

			} );

		}

		return Promise.resolve( sourceURI ).then( function ( sourceURI ) {

			// Load Texture resource.

			loader = options.manager.getHandler( sourceURI ) || loader || parser.textureLoader;

			return new Promise( function ( resolve, reject ) {

				loader.load( resolveURL( sourceURI, options.path ), resolve, undefined, reject );

			} ).then( function ( texture ) {

				// Clean up resources.
				if ( isObjectURL ) URL.revokeObjectURL( sourceURI );

				// Ignore unknown mime types, like DDS files.
				if ( imageDef.mimeType in MIME_TYPE_FORMATS ) texture.format = MIME_TYPE_FORMATS[ imageDef.mimeType ];

				return texture;

			} );

		} );

	};

	/**
	 * Specification: https://github.com/KhronosGroup/glTF/tree/master/specification/2.0#textures
	 * @param {GLTF.definition} textureDef
	 * @return {Promise<THREE.Texture>}
	 */
	GLTFParser.prototype.loadTexture = function ( textureDef ) {

		var json = this.json;

		var pending = this.loadTextureFile( json.images[ textureDef.source ] );

		return pending.then( function ( texture ) {

			texture.flipY = false;

			if ( textureDef.name !== undefined ) texture.name = textureDef.name;

			var samplers = json.samplers || {};
			var sampler = samplers[ textureDef.sampler ] || {};

			texture.magFilter = WEBGL_FILTERS[ sampler.magFilter ] || LinearFilter;
			texture.minFilter = WEBGL_FILTERS[ sampler.minFilter ] || LinearMipmapLinearFilter;
			texture.wrapS = WEBGL_WRAPPINGS[ sampler.wrapS ] || RepeatWrapping;
			texture.wrapT = WEBGL_WRAPPINGS[ sampler.wrapT ] || RepeatWrapping;

			return texture;

		} );

	};

	/**
	 * Asynchronously assigns a texture to the given material parameters.
	 * @param {Object} materialParams
	 * @param {string} mapName
	 * @param {Object} mapDef
	 * @return {Promise}
	 */
	GLTFParser.prototype.assignTexture = function ( materialParams, mapName, mapDef ) {

		var parser = this;

		return this.getDependency( 'texture', mapDef.index ).then( function ( texture ) {

			if ( ! texture.isCompressedTexture ) {

				switch ( mapName ) {

					case 'aoMap':
					case 'emissiveMap':
					case 'metalnessMap':
					case 'normalMap':
					case 'roughnessMap':
						texture.format = RGBFormat;
						break;

				}

			}

			// Materials sample aoMap from UV set 1 and other maps from UV set 0 - this can't be configured
			// However, we will copy UV set 0 to UV set 1 on demand for aoMap
			if ( mapDef.texCoord !== undefined && mapDef.texCoord != 0 && ! ( mapName === 'aoMap' && mapDef.texCoord == 1 ) ) {

				console.warn( 'THREE.GLTFLoader: Custom UV set ' + mapDef.texCoord + ' for texture ' + mapName + ' not yet supported.' );

			}

			return parser._onAfter( 'Map', texture, mapDef );

		} ).then( function ( texture ) {

			materialParams[ mapName ] = texture;

		} );

	};

	/**
	 * Assigns final material to a Mesh, Line, or Points instance. The instance
	 * already has a material (generated from the glTF material options alone)
	 * but reuse of the same glTF material may require multiple threejs materials
	 * to accomodate different primitive types, defines, etc. New materials will
	 * be created if necessary, and reused from a cache.
	 * @param  {Object3D} mesh Mesh, Line, or Points instance.
	 */
	GLTFParser.prototype.assignFinalMaterial = function ( mesh ) {

		var geometry = mesh.geometry;
		var material = mesh.material;
		var extensions = this.plugins.extensions;

		var useVertexTangents = geometry.attributes.tangent !== undefined;
		var useVertexColors = geometry.attributes.color !== undefined;
		var useFlatShading = geometry.attributes.normal === undefined;
		var useSkinning = mesh.isSkinnedMesh === true;
		var useMorphTargets = Object.keys( geometry.morphAttributes ).length > 0;
		var useMorphNormals = useMorphTargets && geometry.morphAttributes.normal !== undefined;

		if ( mesh.isPoints ) {

			var cacheKey = 'PointsMaterial:' + material.uuid;

			var pointsMaterial = this.cache.get( cacheKey );

			if ( ! pointsMaterial ) {

				pointsMaterial = new PointsMaterial();
				Material.prototype.copy.call( pointsMaterial, material );
				pointsMaterial.color.copy( material.color );
				pointsMaterial.map = material.map;
				pointsMaterial.sizeAttenuation = false; // glTF spec says points should be 1px

				this.cache.add( cacheKey, pointsMaterial );

			}

			material = pointsMaterial;

		} else if ( mesh.isLine ) {

			var cacheKey = 'LineBasicMaterial:' + material.uuid;

			var lineMaterial = this.cache.get( cacheKey );

			if ( ! lineMaterial ) {

				lineMaterial = new LineBasicMaterial();
				Material.prototype.copy.call( lineMaterial, material );
				lineMaterial.color.copy( material.color );

				this.cache.add( cacheKey, lineMaterial );

			}

			material = lineMaterial;

		}

		// Clone the material if it will be modified
		if ( useVertexTangents || useVertexColors || useFlatShading || useSkinning || useMorphTargets ) {

			var cacheKey = 'ClonedMaterial:' + material.uuid + ':';

			if ( material.isGLTFMeshStandardSGMaterial ) cacheKey += 'specular-glossiness:';
			if ( useSkinning ) cacheKey += 'skinning:';
			if ( useVertexTangents ) cacheKey += 'vertex-tangents:';
			if ( useVertexColors ) cacheKey += 'vertex-colors:';
			if ( useFlatShading ) cacheKey += 'flat-shading:';
			if ( useMorphTargets ) cacheKey += 'morph-targets:';
			if ( useMorphNormals ) cacheKey += 'morph-normals:';

			var cachedMaterial = this.cache.get( cacheKey );

			if ( ! cachedMaterial ) {

				cachedMaterial = material.clone();

				if ( useSkinning ) cachedMaterial.skinning = true;
				if ( useVertexTangents ) cachedMaterial.vertexTangents = true;
				if ( useVertexColors ) cachedMaterial.vertexColors = VertexColors;
				if ( useFlatShading ) cachedMaterial.flatShading = true;
				if ( useMorphTargets ) cachedMaterial.morphTargets = true;
				if ( useMorphNormals ) cachedMaterial.morphNormals = true;

				this.cache.add( cacheKey, cachedMaterial );

			}

			material = cachedMaterial;

		}

		// workarounds for mesh and geometry

		if ( material.aoMap && geometry.attributes.uv2 === undefined && geometry.attributes.uv !== undefined ) {

			geometry.setAttribute( 'uv2', new BufferAttribute( geometry.attributes.uv.array, 2 ) );

		}

		// https://github.com/mrdoob/three.js/issues/11438#issuecomment-507003995
		if ( material.normalScale && ! useVertexTangents ) {

			material.normalScale.y = - material.normalScale.y;

		}

		mesh.material = material;

	};

	/**
	 * Specification: https://github.com/KhronosGroup/glTF/blob/master/specification/2.0/README.md#materials
	 * @param {number} materialIndex
	 * @param {GLTF.definition} materialDef
	 * @param {Material|null} customMaterial (optional)
	 * @return {Promise<Material>}
	 */
	GLTFParser.prototype.loadMaterial = function ( materialIndex, materialDef, customMaterial ) {

		var parser = this;
		var json = this.json;
		var materialDef = json.materials[ materialIndex ];
		var extensions = this.plugins.extensions;

		var materialPending;

		if ( customMaterial ) {

			materialPending = Promise.resolve( customMaterial );

		} else {

			var materialParams = {};
			var pending = [];

			// Specification:
			// https://github.com/KhronosGroup/glTF/tree/master/specification/2.0#metallic-roughness-material

			var metallicRoughness = materialDef.pbrMetallicRoughness || {};

			materialParams.color = new Color( 1.0, 1.0, 1.0 );
			materialParams.opacity = 1.0;

			if ( Array.isArray( metallicRoughness.baseColorFactor ) ) {

				var array = metallicRoughness.baseColorFactor;

				materialParams.color.fromArray( array );
				materialParams.opacity = array[ 3 ];

			}

			if ( metallicRoughness.baseColorTexture !== undefined ) {

				pending.push( parser.assignTexture( materialParams, 'map', metallicRoughness.baseColorTexture ) );

			}

			materialParams.metalness = metallicRoughness.metallicFactor !== undefined ? metallicRoughness.metallicFactor : 1.0;
			materialParams.roughness = metallicRoughness.roughnessFactor !== undefined ? metallicRoughness.roughnessFactor : 1.0;

			if ( metallicRoughness.metallicRoughnessTexture !== undefined ) {

				pending.push( parser.assignTexture( materialParams, 'metalnessMap', metallicRoughness.metallicRoughnessTexture ) );
				pending.push( parser.assignTexture( materialParams, 'roughnessMap', metallicRoughness.metallicRoughnessTexture ) );

			}

			if ( materialDef.normalTexture !== undefined ) {

				pending.push( parser.assignTexture( materialParams, 'normalMap', materialDef.normalTexture ) );

				materialParams.normalScale = new Vector2( 1, 1 );

				if ( materialDef.normalTexture.scale !== undefined ) {

					materialParams.normalScale.set( materialDef.normalTexture.scale, materialDef.normalTexture.scale );

				}

			}

			if ( materialDef.occlusionTexture !== undefined ) {

				pending.push( parser.assignTexture( materialParams, 'aoMap', materialDef.occlusionTexture ) );

				if ( materialDef.occlusionTexture.strength !== undefined ) {

					materialParams.aoMapIntensity = materialDef.occlusionTexture.strength;

				}

			}

			if ( materialDef.emissiveFactor !== undefined ) {

				materialParams.emissive = new Color().fromArray( materialDef.emissiveFactor );

			}

			if ( materialDef.emissiveTexture !== undefined ) {

				pending.push( parser.assignTexture( materialParams, 'emissiveMap', materialDef.emissiveTexture ) );

			}

			materialPending = Promise.all( pending ).then( function () {

				return new MeshStandardMaterial( materialParams );

			} );

		}

		return materialPending.then( function ( material ) {

			if ( materialDef.name !== undefined ) material.name = materialDef.name;
			if ( materialDef.doubleSided === true ) material.side = DoubleSide;

			var alphaMode = materialDef.alphaMode || ALPHA_MODES.OPAQUE;

			if ( alphaMode === ALPHA_MODES.BLEND ) {

				material.transparent = true;

			} else {

				material.transparent = false;

				if ( alphaMode === ALPHA_MODES.MASK ) {

					material.alphaTest = materialDef.alphaCutoff !== undefined ? materialDef.alphaCutoff : 0.5;

				}

			}

			// baseColorTexture, emissiveTexture, and specularGlossinessTexture use sRGB encoding.
			if ( material.map ) material.map.encoding = sRGBEncoding;
			if ( material.emissiveMap ) material.emissiveMap.encoding = sRGBEncoding;

			assignExtrasToUserData( material, materialDef );

			if ( materialDef.extensions ) addUnknownExtensionsToUserData( extensions, material, materialDef );

			return material;

		} );

	};

	/**
	 * @param {BufferGeometry} geometry
	 * @param {GLTF.Primitive} primitiveDef
	 * @param {GLTFParser} parser
	 */
	function computeBounds( geometry, primitiveDef, parser ) {

		var attributes = primitiveDef.attributes;

		var box = new Box3();

		if ( attributes.POSITION !== undefined ) {

			var accessor = parser.json.accessors[ attributes.POSITION ];

			var min = accessor.min;
			var max = accessor.max;

			// glTF requires 'min' and 'max', but VRM (which extends glTF) currently ignores that requirement.

			if ( min !== undefined && max !== undefined ) {

				box.set(
					new Vector3( min[ 0 ], min[ 1 ], min[ 2 ] ),
					new Vector3( max[ 0 ], max[ 1 ], max[ 2 ] ) );

			} else {

				console.warn( 'THREE.GLTFLoader: Missing min/max properties for accessor POSITION.' );

				return;

			}

		} else {

			return;

		}

		var targets = primitiveDef.targets;

		if ( targets !== undefined ) {

			var vector = new Vector3();

			for ( var i = 0, il = targets.length; i < il; i ++ ) {

				var target = targets[ i ];

				if ( target.POSITION !== undefined ) {

					var accessor = parser.json.accessors[ target.POSITION ];
					var min = accessor.min;
					var max = accessor.max;

					// glTF requires 'min' and 'max', but VRM (which extends glTF) currently ignores that requirement.

					if ( min !== undefined && max !== undefined ) {

						// we need to get max of absolute components because target weight is [-1,1]
						vector.setX( Math.max( Math.abs( min[ 0 ] ), Math.abs( max[ 0 ] ) ) );
						vector.setY( Math.max( Math.abs( min[ 1 ] ), Math.abs( max[ 1 ] ) ) );
						vector.setZ( Math.max( Math.abs( min[ 2 ] ), Math.abs( max[ 2 ] ) ) );

						box.expandByVector( vector );

					} else {

						console.warn( 'THREE.GLTFLoader: Missing min/max properties for accessor POSITION.' );

					}

				}

			}

		}

		geometry.boundingBox = box;

		var sphere = new Sphere();

		box.getCenter( sphere.center );
		sphere.radius = box.min.distanceTo( box.max ) / 2;

		geometry.boundingSphere = sphere;

	}

	/**
	 * @param {BufferGeometry} geometry
	 * @param {GLTF.Primitive} primitiveDef
	 * @param {GLTFParser} parser
	 * @return {Promise<BufferGeometry>}
	 */
	function addPrimitiveAttributes( geometry, primitiveDef, parser ) {

		var attributes = primitiveDef.attributes;

		var pending = [];

		function assignAttributeAccessor( accessorIndex, attributeName ) {

			return parser.getDependency( 'accessor', accessorIndex )
				.then( function ( accessor ) {

					geometry.setAttribute( attributeName, accessor );

				} );

		}

		for ( var gltfAttributeName in attributes ) {

			var threeAttributeName = ATTRIBUTES[ gltfAttributeName ] || gltfAttributeName.toLowerCase();

			// Skip attributes already provided by e.g. Draco extension.
			if ( threeAttributeName in geometry.attributes ) continue;

			pending.push( assignAttributeAccessor( attributes[ gltfAttributeName ], threeAttributeName ) );

		}

		if ( primitiveDef.indices !== undefined && ! geometry.index ) {

			var accessor = parser.getDependency( 'accessor', primitiveDef.indices ).then( function ( accessor ) {

				geometry.setIndex( accessor );

			} );

			pending.push( accessor );

		}

		assignExtrasToUserData( geometry, primitiveDef );

		computeBounds( geometry, primitiveDef, parser );

		return Promise.all( pending ).then( function () {

			return primitiveDef.targets !== undefined
				? addMorphTargets( geometry, primitiveDef.targets, parser )
				: geometry;

		} );

	}

	/**
	 * @param {BufferGeometry} geometry
	 * @param {Number} drawMode
	 * @return {BufferGeometry}
	 */
	function toTrianglesDrawMode( geometry, drawMode ) {

		var index = geometry.getIndex();

		// generate index if not present

		if ( index === null ) {

			var indices = [];

			var position = geometry.getAttribute( 'position' );

			if ( position !== undefined ) {

				for ( var i = 0; i < position.count; i ++ ) {

					indices.push( i );

				}

				geometry.setIndex( indices );
				index = geometry.getIndex();

			} else {

				console.error( 'THREE.GLTFLoader.toTrianglesDrawMode(): Undefined position attribute. Processing not possible.' );
				return geometry;

			}

		}

		//

		var numberOfTriangles = index.count - 2;
		var newIndices = [];

		if ( drawMode === TriangleFanDrawMode ) {

			// gl.TRIANGLE_FAN

			for ( var i = 1; i <= numberOfTriangles; i ++ ) {

				newIndices.push( index.getX( 0 ) );
				newIndices.push( index.getX( i ) );
				newIndices.push( index.getX( i + 1 ) );

			}

		} else {

			// gl.TRIANGLE_STRIP

			for ( var i = 0; i < numberOfTriangles; i ++ ) {

				if ( i % 2 === 0 ) {

					newIndices.push( index.getX( i ) );
					newIndices.push( index.getX( i + 1 ) );
					newIndices.push( index.getX( i + 2 ) );


				} else {

					newIndices.push( index.getX( i + 2 ) );
					newIndices.push( index.getX( i + 1 ) );
					newIndices.push( index.getX( i ) );

				}

			}

		}

		if ( ( newIndices.length / 3 ) !== numberOfTriangles ) {

			console.error( 'THREE.GLTFLoader.toTrianglesDrawMode(): Unable to generate correct amount of triangles.' );

		}

		// build final geometry

		var newGeometry = geometry.clone();
		newGeometry.setIndex( newIndices );

		return newGeometry;

	}

	/**
	 * Specification: https://github.com/KhronosGroup/glTF/blob/master/specification/2.0/README.md#geometry
	 *
	 * Creates BufferGeometries from primitives.
	 *
	 * @param {Array<GLTF.Primitive>} primitives
	 * @return {Promise<Array<BufferGeometry>>}
	 */
	GLTFParser.prototype.loadGeometries = function ( primitives ) {

		var parser = this;
		var extensions = this.plugins.extensions;
		var cache = this.primitiveCache;

		var pending = [];

		for ( var i = 0, il = primitives.length; i < il; i ++ ) {

			var primitive = primitives[ i ];
			var cacheKey = createPrimitiveKey( primitive );

			// See if we've already created this geometry
			var cached = cache[ cacheKey ];

			if ( cached ) {

				// Use the cached geometry if it exists
				pending.push( cached.promise );

			} else {

				var geometryPromise = this._onBefore( 'Geometry', primitive ).then( function ( primitive ) {

					return parser._on( 'Geometry', primitive ).then( function ( geometry ) {

						return geometry || addPrimitiveAttributes( new BufferGeometry(), primitive, parser );

					} ).then( function ( geometry ) {

						return parser._onAfter( 'Geometry', geometry, primitive );

					} );

				} );

				// Cache this geometry
				cache[ cacheKey ] = { primitive: primitive, promise: geometryPromise };

				pending.push( geometryPromise );

			}

		}

		return Promise.all( pending );

	};

	/**
	 * Specification: https://github.com/KhronosGroup/glTF/blob/master/specification/2.0/README.md#meshes
	 * @param {number} meshIndex
	 * @param {GLTF.definition} meshDef
	 * @param {Group|Mesh|SkinnedMesh|null} customMesh (optional)
	 * @return {Promise<Group|Mesh|SkinnedMesh>}
	 */
	GLTFParser.prototype.loadMesh = function ( meshIndex, meshDef, customMesh ) {

		var parser = this;
		var json = this.json;

		var meshDef = json.meshes[ meshIndex ];
		var primitives = meshDef.primitives;

		var pending = [];

		for ( var i = 0, il = primitives.length; i < il; i ++ ) {

			var material = primitives[ i ].material === undefined
				? createDefaultMaterial( this.cache )
				: this.getDependency( 'material', primitives[ i ].material );

			pending.push( material );

		}

		pending.push( parser.loadGeometries( primitives ) );

		return Promise.all( pending ).then( function ( results ) {

			var materials = results.slice( 0, results.length - 1 );
			var geometries = results[ results.length - 1 ];

			var meshes = [];

			for ( var i = 0, il = geometries.length; i < il; i ++ ) {

				var geometry = geometries[ i ];
				var primitive = primitives[ i ];

				// 1. create Mesh

				var mesh;

				var material = materials[ i ];

				if ( primitive.mode === WEBGL_CONSTANTS.TRIANGLES ||
					primitive.mode === WEBGL_CONSTANTS.TRIANGLE_STRIP ||
					primitive.mode === WEBGL_CONSTANTS.TRIANGLE_FAN ||
					primitive.mode === undefined ) {

					// .isSkinnedMesh isn't in glTF spec. See .markDefs()
					mesh = meshDef.isSkinnedMesh === true
						? new SkinnedMesh( geometry, material )
						: new Mesh( geometry, material );

					if ( mesh.isSkinnedMesh === true && ! mesh.geometry.attributes.skinWeight.normalized ) {

						// we normalize floating point skin weight array to fix malformed assets (see #15319)
						// it's important to skip this for non-float32 data since normalizeSkinWeights assumes non-normalized inputs
						mesh.normalizeSkinWeights();

					}

					if ( primitive.mode === WEBGL_CONSTANTS.TRIANGLE_STRIP ) {

						mesh.geometry = toTrianglesDrawMode( mesh.geometry, TriangleStripDrawMode );

					} else if ( primitive.mode === WEBGL_CONSTANTS.TRIANGLE_FAN ) {

						mesh.geometry = toTrianglesDrawMode( mesh.geometry, TriangleFanDrawMode );

					}

				} else if ( primitive.mode === WEBGL_CONSTANTS.LINES ) {

					mesh = new LineSegments( geometry, material );

				} else if ( primitive.mode === WEBGL_CONSTANTS.LINE_STRIP ) {

					mesh = new Line( geometry, material );

				} else if ( primitive.mode === WEBGL_CONSTANTS.LINE_LOOP ) {

					mesh = new LineLoop( geometry, material );

				} else if ( primitive.mode === WEBGL_CONSTANTS.POINTS ) {

					mesh = new Points( geometry, material );

				} else {

					throw new Error( 'THREE.GLTFLoader: Primitive mode unsupported: ' + primitive.mode );

				}

				if ( Object.keys( mesh.geometry.morphAttributes ).length > 0 ) {

					updateMorphTargets( mesh, meshDef );

				}

				mesh.name = meshDef.name || ( 'mesh_' + meshIndex );

				if ( geometries.length > 1 ) mesh.name += '_' + i;

				assignExtrasToUserData( mesh, meshDef );

				parser.assignFinalMaterial( mesh );

				meshes.push( mesh );

			}

			if ( meshes.length === 1 ) {

				return meshes[ 0 ];

			}

			var group = new Group();

			for ( var i = 0, il = meshes.length; i < il; i ++ ) {

				group.add( meshes[ i ] );

			}

			return group;

		} );

	};

	/**
	 * Specification: https://github.com/KhronosGroup/glTF/tree/master/specification/2.0#cameras
	 * @param {GLTF.definition} cameraDef
	 * @return {Promise<THREE.Camera>}
	 */
	GLTFParser.prototype.loadCamera = function ( cameraDef ) {

		var camera;
		var params = cameraDef[ cameraDef.type ];

		if ( ! params ) {

			console.warn( 'THREE.GLTFLoader: Missing camera parameters.' );
			return;

		}

		if ( cameraDef.type === 'perspective' ) {

			camera = new PerspectiveCamera( MathUtils.radToDeg( params.yfov ), params.aspectRatio || 1, params.znear || 1, params.zfar || 2e6 );

		} else if ( cameraDef.type === 'orthographic' ) {

			camera = new OrthographicCamera( params.xmag / - 2, params.xmag / 2, params.ymag / 2, params.ymag / - 2, params.znear, params.zfar );

		}

		if ( cameraDef.name !== undefined ) camera.name = cameraDef.name;

		assignExtrasToUserData( camera, cameraDef );

		return Promise.resolve( camera );

	};

	/**
	 * Specification: https://github.com/KhronosGroup/glTF/tree/master/specification/2.0#skins
	 * @param {GLTF.definition} skinDef
	 * @return {Promise<Object>}
	 */
	GLTFParser.prototype.loadSkin = function ( skinDef ) {

		var skinEntry = { joints: skinDef.joints };

		if ( skinDef.inverseBindMatrices === undefined ) {

			return Promise.resolve( skinEntry );

		}

		return this.getDependency( 'accessor', skinDef.inverseBindMatrices ).then( function ( accessor ) {

			skinEntry.inverseBindMatrices = accessor;

			return skinEntry;

		} );

	};

	/**
	 * Specification: https://github.com/KhronosGroup/glTF/tree/master/specification/2.0#animations
	 * @param {number} animationIndex
	 * @param {GLTF.definition} animationDef
	 * @return {Promise<AnimationClip>}
	 */
	GLTFParser.prototype.loadAnimation = function ( animationIndex, animationDef ) {

		var json = this.json;

		var pendingNodes = [];
		var pendingInputAccessors = [];
		var pendingOutputAccessors = [];
		var pendingSamplers = [];
		var pendingTargets = [];

		for ( var i = 0, il = animationDef.channels.length; i < il; i ++ ) {

			var channel = animationDef.channels[ i ];
			var sampler = animationDef.samplers[ channel.sampler ];
			var target = channel.target;
			var name = target.node !== undefined ? target.node : target.id; // NOTE: target.id is deprecated.
			var input = animationDef.parameters !== undefined ? animationDef.parameters[ sampler.input ] : sampler.input;
			var output = animationDef.parameters !== undefined ? animationDef.parameters[ sampler.output ] : sampler.output;

			pendingNodes.push( this.getDependency( 'node', name ) );
			pendingInputAccessors.push( this.getDependency( 'accessor', input ) );
			pendingOutputAccessors.push( this.getDependency( 'accessor', output ) );
			pendingSamplers.push( sampler );
			pendingTargets.push( target );

		}

		return Promise.all( [

			Promise.all( pendingNodes ),
			Promise.all( pendingInputAccessors ),
			Promise.all( pendingOutputAccessors ),
			Promise.all( pendingSamplers ),
			Promise.all( pendingTargets )

		] ).then( function ( dependencies ) {

			var nodes = dependencies[ 0 ];
			var inputAccessors = dependencies[ 1 ];
			var outputAccessors = dependencies[ 2 ];
			var samplers = dependencies[ 3 ];
			var targets = dependencies[ 4 ];

			var tracks = [];

			for ( var i = 0, il = nodes.length; i < il; i ++ ) {

				var node = nodes[ i ];
				var inputAccessor = inputAccessors[ i ];
				var outputAccessor = outputAccessors[ i ];
				var sampler = samplers[ i ];
				var target = targets[ i ];

				if ( node === undefined ) continue;

				node.updateMatrix();
				node.matrixAutoUpdate = true;

				var TypedKeyframeTrack;

				switch ( PATH_PROPERTIES[ target.path ] ) {

					case PATH_PROPERTIES.weights:

						TypedKeyframeTrack = NumberKeyframeTrack;
						break;

					case PATH_PROPERTIES.rotation:

						TypedKeyframeTrack = QuaternionKeyframeTrack;
						break;

					case PATH_PROPERTIES.position:
					case PATH_PROPERTIES.scale:
					default:

						TypedKeyframeTrack = VectorKeyframeTrack;
						break;

				}

				var targetName = node.name ? node.name : node.uuid;

				var interpolation = sampler.interpolation !== undefined ? INTERPOLATION[ sampler.interpolation ] : InterpolateLinear;

				var targetNames = [];

				if ( PATH_PROPERTIES[ target.path ] === PATH_PROPERTIES.weights ) {

					// Node may be a Group (glTF mesh with several primitives) or a Mesh.
					node.traverse( function ( object ) {

						if ( object.isMesh === true && object.morphTargetInfluences ) {

							targetNames.push( object.name ? object.name : object.uuid );

						}

					} );

				} else {

					targetNames.push( targetName );

				}

				var outputArray = outputAccessor.array;

				if ( outputAccessor.normalized ) {

					var scale;

					if ( outputArray.constructor === Int8Array ) {

						scale = 1 / 127;

					} else if ( outputArray.constructor === Uint8Array ) {

						scale = 1 / 255;

					} else if ( outputArray.constructor == Int16Array ) {

						scale = 1 / 32767;

					} else if ( outputArray.constructor === Uint16Array ) {

						scale = 1 / 65535;

					} else {

						throw new Error( 'THREE.GLTFLoader: Unsupported output accessor component type.' );

					}

					var scaled = new Float32Array( outputArray.length );

					for ( var j = 0, jl = outputArray.length; j < jl; j ++ ) {

						scaled[ j ] = outputArray[ j ] * scale;

					}

					outputArray = scaled;

				}

				for ( var j = 0, jl = targetNames.length; j < jl; j ++ ) {

					var track = new TypedKeyframeTrack(
						targetNames[ j ] + '.' + PATH_PROPERTIES[ target.path ],
						inputAccessor.array,
						outputArray,
						interpolation
					);

					// Override interpolation with custom factory method.
					if ( sampler.interpolation === 'CUBICSPLINE' ) {

						track.createInterpolant = function InterpolantFactoryMethodGLTFCubicSpline( result ) {

							// A CUBICSPLINE keyframe in glTF has three output values for each input value,
							// representing inTangent, splineVertex, and outTangent. As a result, track.getValueSize()
							// must be divided by three to get the interpolant's sampleSize argument.

							return new GLTFCubicSplineInterpolant( this.times, this.values, this.getValueSize() / 3, result );

						};

						// Mark as CUBICSPLINE. `track.getInterpolation()` doesn't support custom interpolants.
						track.createInterpolant.isInterpolantFactoryMethodGLTFCubicSpline = true;

					}

					tracks.push( track );

				}

			}

			var name = animationDef.name !== undefined ? animationDef.name : 'animation_' + animationIndex;

			return new AnimationClip( name, undefined, tracks );

		} );

	};

	/**
	 * Specification: https://github.com/KhronosGroup/glTF/tree/master/specification/2.0#nodes-and-hierarchy
	 * @param {number} nodeIndex
	 * @param {GLTF.definition} nodeDef
	 * @param {Object3D|null} customNode (optional)
	 * @return {Promise<Object3D>}
	 */
	GLTFParser.prototype.loadNode = function ( nodeIndex, nodeDef, customNode ) {

		var json = this.json;
		var extensions = this.plugins.extensions;
		var parser = this;

		var meshReferences = json.meshReferences;
		var meshUses = json.meshUses;

		return ( function () {

			var pending = [];

			if ( customNode ) {

				pending.push( Promise.resolve( customNode ) );

			}

			if ( nodeDef.mesh !== undefined ) {

				pending.push( parser.getDependency( 'mesh', nodeDef.mesh ).then( function ( mesh ) {

					var node;

					if ( meshReferences[ nodeDef.mesh ] > 1 ) {

						var instanceNum = meshUses[ nodeDef.mesh ] ++;

						node = mesh.clone();
						node.name += '_instance_' + instanceNum;

					} else {

						node = mesh;

					}

					// if weights are provided on the node, override weights on the mesh.
					if ( nodeDef.weights !== undefined ) {

						node.traverse( function ( o ) {

							if ( ! o.isMesh ) return;

							for ( var i = 0, il = nodeDef.weights.length; i < il; i ++ ) {

								o.morphTargetInfluences[ i ] = nodeDef.weights[ i ];

							}

						} );

					}

					return node;

				} ) );

			}

			if ( nodeDef.camera !== undefined ) {

				pending.push( parser.getDependency( 'camera', nodeDef.camera ) );

			}

			return Promise.all( pending );

		}() ).then( function ( objects ) {

			var node;

			// .isBone isn't in glTF spec. See .markDefs
			if ( nodeDef.isBone === true ) {

				node = new Bone();

			} else if ( objects.length > 1 ) {

				node = new Group();

			} else if ( objects.length === 1 ) {

				node = objects[ 0 ];

			} else {

				node = new Object3D();

			}

			if ( node !== objects[ 0 ] ) {

				for ( var i = 0, il = objects.length; i < il; i ++ ) {

					node.add( objects[ i ] );

				}

			}

			if ( nodeDef.name !== undefined ) {

				node.userData.name = nodeDef.name;
				node.name = PropertyBinding.sanitizeNodeName( nodeDef.name );

			}

			assignExtrasToUserData( node, nodeDef );

			if ( nodeDef.extensions ) addUnknownExtensionsToUserData( extensions, node, nodeDef );

			if ( nodeDef.matrix !== undefined ) {

				var matrix = new Matrix4();
				matrix.fromArray( nodeDef.matrix );
				node.applyMatrix( matrix );

			} else {

				if ( nodeDef.translation !== undefined ) {

					node.position.fromArray( nodeDef.translation );

				}

				if ( nodeDef.rotation !== undefined ) {

					node.quaternion.fromArray( nodeDef.rotation );

				}

				if ( nodeDef.scale !== undefined ) {

					node.scale.fromArray( nodeDef.scale );

				}

			}

			return node;

		} );

	};

	/**
	 * Specification: https://github.com/KhronosGroup/glTF/tree/master/specification/2.0#scenes
	 * @param {number} sceneIndex
	 * @param {GLTF.definition} sceneDef
	 * @param {Scene|null} customScene (optional)
	 * @return {Promise<Scene>}
	 */
	GLTFParser.prototype.loadScene = function () {

		// scene node hierachy builder

		function buildNodeHierachy( nodeId, parentObject, json, parser ) {

			var nodeDef = json.nodes[ nodeId ];

			return parser.getDependency( 'node', nodeId ).then( function ( node ) {

				if ( nodeDef.skin === undefined ) return node;

				// build skeleton here as well

				var skinEntry;

				return parser.getDependency( 'skin', nodeDef.skin ).then( function ( skin ) {

					skinEntry = skin;

					var pendingJoints = [];

					for ( var i = 0, il = skinEntry.joints.length; i < il; i ++ ) {

						pendingJoints.push( parser.getDependency( 'node', skinEntry.joints[ i ] ) );

					}

					return Promise.all( pendingJoints );

				} ).then( function ( jointNodes ) {

					node.traverse( function ( mesh ) {

						if ( ! mesh.isMesh ) return;

						var bones = [];
						var boneInverses = [];

						for ( var j = 0, jl = jointNodes.length; j < jl; j ++ ) {

							var jointNode = jointNodes[ j ];

							if ( jointNode ) {

								bones.push( jointNode );

								var mat = new Matrix4();

								if ( skinEntry.inverseBindMatrices !== undefined ) {

									mat.fromArray( skinEntry.inverseBindMatrices.array, j * 16 );

								}

								boneInverses.push( mat );

							} else {

								console.warn( 'THREE.GLTFLoader: Joint "%s" could not be found.', skinEntry.joints[ j ] );

							}

						}

						mesh.bind( new Skeleton( bones, boneInverses ), mesh.matrixWorld );

					} );

					return node;

				} );

			} ).then( function ( node ) {

				// build node hierachy

				parentObject.add( node );

				var pending = [];

				if ( nodeDef.children ) {

					var children = nodeDef.children;

					for ( var i = 0, il = children.length; i < il; i ++ ) {

						var child = children[ i ];
						pending.push( buildNodeHierachy( child, node, json, parser ) );

					}

				}

				return Promise.all( pending );

			} );

		}

		return function loadScene( sceneIndex, sceneDef, customScene ) {

			var json = this.json;
			var extensions = this.plugins.extensions;
			var parser = this;

			var scenePending;

			if ( customScene ) {

				scenePending = Promise.resolve( customScene );

			} else {

				var scene = new Scene();

				var nodeIds = sceneDef.nodes || [];

				var pending = [];

				for ( var i = 0, il = nodeIds.length; i < il; i ++ ) {

					pending.push( buildNodeHierachy( nodeIds[ i ], scene, json, parser ) );

				}

				scenePending = Promise.all( pending ).then( function () {

					return scene;

				} );

			}

			return scenePending.then( function ( scene ) {

				if ( sceneDef.name !== undefined ) scene.name = sceneDef.name;

				assignExtrasToUserData( scene, sceneDef );

				if ( sceneDef.extensions ) addUnknownExtensionsToUserData( extensions, scene, sceneDef );

				return scene;

			} );

		};

	}();

	return GLTFLoader;

} )();

export { GLTFLoader };
