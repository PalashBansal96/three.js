// threejs.org/license
(function (global, factory) {
	typeof exports === 'object' && typeof module !== 'undefined' ? factory(exports) :
	typeof define === 'function' && define.amd ? define(['exports'], factory) :
	(global = typeof globalThis !== 'undefined' ? globalThis : global || self, factory(global.THREE_VEC = {}));
}(this, (function (exports) { 'use strict';

	var _lut = [];

	for (var i = 0; i < 256; i++) {
		_lut[i] = (i < 16 ? '0' : '') + i.toString(16);
	}

	var _seed = 1234567;
	var MathUtils = {
		DEG2RAD: Math.PI / 180,
		RAD2DEG: 180 / Math.PI,
		generateUUID: function generateUUID() {
			// http://stackoverflow.com/questions/105034/how-to-create-a-guid-uuid-in-javascript/21963136#21963136
			var d0 = Math.random() * 0xffffffff | 0;
			var d1 = Math.random() * 0xffffffff | 0;
			var d2 = Math.random() * 0xffffffff | 0;
			var d3 = Math.random() * 0xffffffff | 0;
			var uuid = _lut[d0 & 0xff] + _lut[d0 >> 8 & 0xff] + _lut[d0 >> 16 & 0xff] + _lut[d0 >> 24 & 0xff] + '-' + _lut[d1 & 0xff] + _lut[d1 >> 8 & 0xff] + '-' + _lut[d1 >> 16 & 0x0f | 0x40] + _lut[d1 >> 24 & 0xff] + '-' + _lut[d2 & 0x3f | 0x80] + _lut[d2 >> 8 & 0xff] + '-' + _lut[d2 >> 16 & 0xff] + _lut[d2 >> 24 & 0xff] + _lut[d3 & 0xff] + _lut[d3 >> 8 & 0xff] + _lut[d3 >> 16 & 0xff] + _lut[d3 >> 24 & 0xff]; // .toUpperCase() here flattens concatenated strings to save heap memory space.

			return uuid.toUpperCase();
		},
		clamp: function clamp(value, min, max) {
			return Math.max(min, Math.min(max, value));
		},
		// compute euclidian modulo of m % n
		// https://en.wikipedia.org/wiki/Modulo_operation
		euclideanModulo: function euclideanModulo(n, m) {
			return (n % m + m) % m;
		},
		// Linear mapping from range <a1, a2> to range <b1, b2>
		mapLinear: function mapLinear(x, a1, a2, b1, b2) {
			return b1 + (x - a1) * (b2 - b1) / (a2 - a1);
		},
		// https://en.wikipedia.org/wiki/Linear_interpolation
		lerp: function lerp(x, y, t) {
			return (1 - t) * x + t * y;
		},
		// http://en.wikipedia.org/wiki/Smoothstep
		smoothstep: function smoothstep(x, min, max) {
			if (x <= min) return 0;
			if (x >= max) return 1;
			x = (x - min) / (max - min);
			return x * x * (3 - 2 * x);
		},
		smootherstep: function smootherstep(x, min, max) {
			if (x <= min) return 0;
			if (x >= max) return 1;
			x = (x - min) / (max - min);
			return x * x * x * (x * (x * 6 - 15) + 10);
		},
		// Random integer from <low, high> interval
		randInt: function randInt(low, high) {
			return low + Math.floor(Math.random() * (high - low + 1));
		},
		// Random float from <low, high> interval
		randFloat: function randFloat(low, high) {
			return low + Math.random() * (high - low);
		},
		// Random float from <-range/2, range/2> interval
		randFloatSpread: function randFloatSpread(range) {
			return range * (0.5 - Math.random());
		},
		// Deterministic pseudo-random float in the interval [ 0, 1 ]
		seededRandom: function seededRandom(s) {
			if (s !== undefined) _seed = s % 2147483647; // Park-Miller algorithm

			_seed = _seed * 16807 % 2147483647;
			return (_seed - 1) / 2147483646;
		},
		degToRad: function degToRad(degrees) {
			return degrees * MathUtils.DEG2RAD;
		},
		radToDeg: function radToDeg(radians) {
			return radians * MathUtils.RAD2DEG;
		},
		isPowerOfTwo: function isPowerOfTwo(value) {
			return (value & value - 1) === 0 && value !== 0;
		},
		ceilPowerOfTwo: function ceilPowerOfTwo(value) {
			return Math.pow(2, Math.ceil(Math.log(value) / Math.LN2));
		},
		floorPowerOfTwo: function floorPowerOfTwo(value) {
			return Math.pow(2, Math.floor(Math.log(value) / Math.LN2));
		},
		setQuaternionFromProperEuler: function setQuaternionFromProperEuler(q, a, b, c, order) {
			// Intrinsic Proper Euler Angles - see https://en.wikipedia.org/wiki/Euler_angles
			// rotations are applied to the axes in the order specified by 'order'
			// rotation by angle 'a' is applied first, then by angle 'b', then by angle 'c'
			// angles are in radians
			var cos = Math.cos;
			var sin = Math.sin;
			var c2 = cos(b / 2);
			var s2 = sin(b / 2);
			var c13 = cos((a + c) / 2);
			var s13 = sin((a + c) / 2);
			var c1_3 = cos((a - c) / 2);
			var s1_3 = sin((a - c) / 2);
			var c3_1 = cos((c - a) / 2);
			var s3_1 = sin((c - a) / 2);

			switch (order) {
				case 'XYX':
					q.set(c2 * s13, s2 * c1_3, s2 * s1_3, c2 * c13);
					break;

				case 'YZY':
					q.set(s2 * s1_3, c2 * s13, s2 * c1_3, c2 * c13);
					break;

				case 'ZXZ':
					q.set(s2 * c1_3, s2 * s1_3, c2 * s13, c2 * c13);
					break;

				case 'XZX':
					q.set(c2 * s13, s2 * s3_1, s2 * c3_1, c2 * c13);
					break;

				case 'YXY':
					q.set(s2 * c3_1, c2 * s13, s2 * s3_1, c2 * c13);
					break;

				case 'ZYZ':
					q.set(s2 * s3_1, s2 * c3_1, c2 * s13, c2 * c13);
					break;

				default:
					console.warn('THREE.MathUtils: .setQuaternionFromProperEuler() encountered an unknown order: ' + order);
			}
		}
	};

	function _defineProperties(target, props) {
		for (var i = 0; i < props.length; i++) {
			var descriptor = props[i];
			descriptor.enumerable = descriptor.enumerable || false;
			descriptor.configurable = true;
			if ("value" in descriptor) descriptor.writable = true;
			Object.defineProperty(target, descriptor.key, descriptor);
		}
	}

	function _createClass(Constructor, protoProps, staticProps) {
		if (protoProps) _defineProperties(Constructor.prototype, protoProps);
		if (staticProps) _defineProperties(Constructor, staticProps);
		return Constructor;
	}

	var Vector2 = /*#__PURE__*/function () {
		function Vector2(x, y) {
			if (x === void 0) {
				x = 0;
			}

			if (y === void 0) {
				y = 0;
			}

			Object.defineProperty(this, 'isVector2', {
				value: true
			});
			this.x = x;
			this.y = y;
		}

		var _proto = Vector2.prototype;

		_proto.set = function set(x, y) {
			this.x = x;
			this.y = y;
			return this;
		};

		_proto.setScalar = function setScalar(scalar) {
			this.x = scalar;
			this.y = scalar;
			return this;
		};

		_proto.setX = function setX(x) {
			this.x = x;
			return this;
		};

		_proto.setY = function setY(y) {
			this.y = y;
			return this;
		};

		_proto.setComponent = function setComponent(index, value) {
			switch (index) {
				case 0:
					this.x = value;
					break;

				case 1:
					this.y = value;
					break;

				default:
					throw new Error('index is out of range: ' + index);
			}

			return this;
		};

		_proto.getComponent = function getComponent(index) {
			switch (index) {
				case 0:
					return this.x;

				case 1:
					return this.y;

				default:
					throw new Error('index is out of range: ' + index);
			}
		};

		_proto.clone = function clone() {
			return new this.constructor(this.x, this.y);
		};

		_proto.copy = function copy(v) {
			this.x = v.x;
			this.y = v.y;
			return this;
		};

		_proto.add = function add(v, w) {
			if (w !== undefined) {
				console.warn('THREE.Vector2: .add() now only accepts one argument. Use .addVectors( a, b ) instead.');
				return this.addVectors(v, w);
			}

			this.x += v.x;
			this.y += v.y;
			return this;
		};

		_proto.addScalar = function addScalar(s) {
			this.x += s;
			this.y += s;
			return this;
		};

		_proto.addVectors = function addVectors(a, b) {
			this.x = a.x + b.x;
			this.y = a.y + b.y;
			return this;
		};

		_proto.addScaledVector = function addScaledVector(v, s) {
			this.x += v.x * s;
			this.y += v.y * s;
			return this;
		};

		_proto.sub = function sub(v, w) {
			if (w !== undefined) {
				console.warn('THREE.Vector2: .sub() now only accepts one argument. Use .subVectors( a, b ) instead.');
				return this.subVectors(v, w);
			}

			this.x -= v.x;
			this.y -= v.y;
			return this;
		};

		_proto.subScalar = function subScalar(s) {
			this.x -= s;
			this.y -= s;
			return this;
		};

		_proto.subVectors = function subVectors(a, b) {
			this.x = a.x - b.x;
			this.y = a.y - b.y;
			return this;
		};

		_proto.multiply = function multiply(v) {
			this.x *= v.x;
			this.y *= v.y;
			return this;
		};

		_proto.multiplyScalar = function multiplyScalar(scalar) {
			this.x *= scalar;
			this.y *= scalar;
			return this;
		};

		_proto.divide = function divide(v) {
			this.x /= v.x;
			this.y /= v.y;
			return this;
		};

		_proto.divideScalar = function divideScalar(scalar) {
			return this.multiplyScalar(1 / scalar);
		};

		_proto.applyMatrix3 = function applyMatrix3(m) {
			var x = this.x,
					y = this.y;
			var e = m.elements;
			this.x = e[0] * x + e[3] * y + e[6];
			this.y = e[1] * x + e[4] * y + e[7];
			return this;
		};

		_proto.min = function min(v) {
			this.x = Math.min(this.x, v.x);
			this.y = Math.min(this.y, v.y);
			return this;
		};

		_proto.max = function max(v) {
			this.x = Math.max(this.x, v.x);
			this.y = Math.max(this.y, v.y);
			return this;
		};

		_proto.clamp = function clamp(min, max) {
			// assumes min < max, componentwise
			this.x = Math.max(min.x, Math.min(max.x, this.x));
			this.y = Math.max(min.y, Math.min(max.y, this.y));
			return this;
		};

		_proto.clampScalar = function clampScalar(minVal, maxVal) {
			this.x = Math.max(minVal, Math.min(maxVal, this.x));
			this.y = Math.max(minVal, Math.min(maxVal, this.y));
			return this;
		};

		_proto.clampLength = function clampLength(min, max) {
			var length = this.length();
			return this.divideScalar(length || 1).multiplyScalar(Math.max(min, Math.min(max, length)));
		};

		_proto.floor = function floor() {
			this.x = Math.floor(this.x);
			this.y = Math.floor(this.y);
			return this;
		};

		_proto.ceil = function ceil() {
			this.x = Math.ceil(this.x);
			this.y = Math.ceil(this.y);
			return this;
		};

		_proto.round = function round() {
			this.x = Math.round(this.x);
			this.y = Math.round(this.y);
			return this;
		};

		_proto.roundToZero = function roundToZero() {
			this.x = this.x < 0 ? Math.ceil(this.x) : Math.floor(this.x);
			this.y = this.y < 0 ? Math.ceil(this.y) : Math.floor(this.y);
			return this;
		};

		_proto.negate = function negate() {
			this.x = -this.x;
			this.y = -this.y;
			return this;
		};

		_proto.dot = function dot(v) {
			return this.x * v.x + this.y * v.y;
		};

		_proto.cross = function cross(v) {
			return this.x * v.y - this.y * v.x;
		};

		_proto.lengthSq = function lengthSq() {
			return this.x * this.x + this.y * this.y;
		};

		_proto.length = function length() {
			return Math.sqrt(this.x * this.x + this.y * this.y);
		};

		_proto.manhattanLength = function manhattanLength() {
			return Math.abs(this.x) + Math.abs(this.y);
		};

		_proto.normalize = function normalize() {
			return this.divideScalar(this.length() || 1);
		};

		_proto.angle = function angle() {
			// computes the angle in radians with respect to the positive x-axis
			var angle = Math.atan2(-this.y, -this.x) + Math.PI;
			return angle;
		};

		_proto.distanceTo = function distanceTo(v) {
			return Math.sqrt(this.distanceToSquared(v));
		};

		_proto.distanceToSquared = function distanceToSquared(v) {
			var dx = this.x - v.x,
					dy = this.y - v.y;
			return dx * dx + dy * dy;
		};

		_proto.manhattanDistanceTo = function manhattanDistanceTo(v) {
			return Math.abs(this.x - v.x) + Math.abs(this.y - v.y);
		};

		_proto.setLength = function setLength(length) {
			return this.normalize().multiplyScalar(length);
		};

		_proto.lerp = function lerp(v, alpha) {
			this.x += (v.x - this.x) * alpha;
			this.y += (v.y - this.y) * alpha;
			return this;
		};

		_proto.lerpVectors = function lerpVectors(v1, v2, alpha) {
			this.x = v1.x + (v2.x - v1.x) * alpha;
			this.y = v1.y + (v2.y - v1.y) * alpha;
			return this;
		};

		_proto.equals = function equals(v) {
			return v.x === this.x && v.y === this.y;
		};

		_proto.fromArray = function fromArray(array, offset) {
			if (offset === undefined) offset = 0;
			this.x = array[offset];
			this.y = array[offset + 1];
			return this;
		};

		_proto.toArray = function toArray(array, offset) {
			if (array === undefined) array = [];
			if (offset === undefined) offset = 0;
			array[offset] = this.x;
			array[offset + 1] = this.y;
			return array;
		};

		_proto.fromBufferAttribute = function fromBufferAttribute(attribute, index, offset) {
			if (offset !== undefined) {
				console.warn('THREE.Vector2: offset has been removed from .fromBufferAttribute().');
			}

			this.x = attribute.getX(index);
			this.y = attribute.getY(index);
			return this;
		};

		_proto.rotateAround = function rotateAround(center, angle) {
			var c = Math.cos(angle),
					s = Math.sin(angle);
			var x = this.x - center.x;
			var y = this.y - center.y;
			this.x = x * c - y * s + center.x;
			this.y = x * s + y * c + center.y;
			return this;
		};

		_proto.random = function random() {
			this.x = Math.random();
			this.y = Math.random();
			return this;
		};

		_createClass(Vector2, [{
			key: "width",
			get: function get() {
				return this.x;
			},
			set: function set(value) {
				this.x = value;
			}
		}, {
			key: "height",
			get: function get() {
				return this.y;
			},
			set: function set(value) {
				this.y = value;
			}
		}]);

		return Vector2;
	}();

	var Quaternion = /*#__PURE__*/function () {
		function Quaternion(x, y, z, w) {
			if (x === void 0) {
				x = 0;
			}

			if (y === void 0) {
				y = 0;
			}

			if (z === void 0) {
				z = 0;
			}

			if (w === void 0) {
				w = 1;
			}

			Object.defineProperty(this, 'isQuaternion', {
				value: true
			});
			this._x = x;
			this._y = y;
			this._z = z;
			this._w = w;
		}

		Quaternion.slerp = function slerp(qa, qb, qm, t) {
			return qm.copy(qa).slerp(qb, t);
		};

		Quaternion.slerpFlat = function slerpFlat(dst, dstOffset, src0, srcOffset0, src1, srcOffset1, t) {
			// fuzz-free, array-based Quaternion SLERP operation
			var x0 = src0[srcOffset0 + 0],
					y0 = src0[srcOffset0 + 1],
					z0 = src0[srcOffset0 + 2],
					w0 = src0[srcOffset0 + 3];
			var x1 = src1[srcOffset1 + 0],
					y1 = src1[srcOffset1 + 1],
					z1 = src1[srcOffset1 + 2],
					w1 = src1[srcOffset1 + 3];

			if (w0 !== w1 || x0 !== x1 || y0 !== y1 || z0 !== z1) {
				var s = 1 - t;
				var cos = x0 * x1 + y0 * y1 + z0 * z1 + w0 * w1,
						dir = cos >= 0 ? 1 : -1,
						sqrSin = 1 - cos * cos; // Skip the Slerp for tiny steps to avoid numeric problems:

				if (sqrSin > Number.EPSILON) {
					var sin = Math.sqrt(sqrSin),
							len = Math.atan2(sin, cos * dir);
					s = Math.sin(s * len) / sin;
					t = Math.sin(t * len) / sin;
				}

				var tDir = t * dir;
				x0 = x0 * s + x1 * tDir;
				y0 = y0 * s + y1 * tDir;
				z0 = z0 * s + z1 * tDir;
				w0 = w0 * s + w1 * tDir; // Normalize in case we just did a lerp:

				if (s === 1 - t) {
					var f = 1 / Math.sqrt(x0 * x0 + y0 * y0 + z0 * z0 + w0 * w0);
					x0 *= f;
					y0 *= f;
					z0 *= f;
					w0 *= f;
				}
			}

			dst[dstOffset] = x0;
			dst[dstOffset + 1] = y0;
			dst[dstOffset + 2] = z0;
			dst[dstOffset + 3] = w0;
		};

		Quaternion.multiplyQuaternionsFlat = function multiplyQuaternionsFlat(dst, dstOffset, src0, srcOffset0, src1, srcOffset1) {
			var x0 = src0[srcOffset0];
			var y0 = src0[srcOffset0 + 1];
			var z0 = src0[srcOffset0 + 2];
			var w0 = src0[srcOffset0 + 3];
			var x1 = src1[srcOffset1];
			var y1 = src1[srcOffset1 + 1];
			var z1 = src1[srcOffset1 + 2];
			var w1 = src1[srcOffset1 + 3];
			dst[dstOffset] = x0 * w1 + w0 * x1 + y0 * z1 - z0 * y1;
			dst[dstOffset + 1] = y0 * w1 + w0 * y1 + z0 * x1 - x0 * z1;
			dst[dstOffset + 2] = z0 * w1 + w0 * z1 + x0 * y1 - y0 * x1;
			dst[dstOffset + 3] = w0 * w1 - x0 * x1 - y0 * y1 - z0 * z1;
			return dst;
		};

		var _proto = Quaternion.prototype;

		_proto.set = function set(x, y, z, w) {
			this._x = x;
			this._y = y;
			this._z = z;
			this._w = w;

			this._onChangeCallback();

			return this;
		};

		_proto.clone = function clone() {
			return new this.constructor(this._x, this._y, this._z, this._w);
		};

		_proto.copy = function copy(quaternion) {
			this._x = quaternion.x;
			this._y = quaternion.y;
			this._z = quaternion.z;
			this._w = quaternion.w;

			this._onChangeCallback();

			return this;
		};

		_proto.setFromEuler = function setFromEuler(euler, update) {
			if (!(euler && euler.isEuler)) {
				throw new Error('THREE.Quaternion: .setFromEuler() now expects an Euler rotation rather than a Vector3 and order.');
			}

			var x = euler._x,
					y = euler._y,
					z = euler._z,
					order = euler._order; // http://www.mathworks.com/matlabcentral/fileexchange/
			// 	20696-function-to-convert-between-dcm-euler-angles-quaternions-and-euler-vectors/
			//	content/SpinCalc.m

			var cos = Math.cos;
			var sin = Math.sin;
			var c1 = cos(x / 2);
			var c2 = cos(y / 2);
			var c3 = cos(z / 2);
			var s1 = sin(x / 2);
			var s2 = sin(y / 2);
			var s3 = sin(z / 2);

			switch (order) {
				case 'XYZ':
					this._x = s1 * c2 * c3 + c1 * s2 * s3;
					this._y = c1 * s2 * c3 - s1 * c2 * s3;
					this._z = c1 * c2 * s3 + s1 * s2 * c3;
					this._w = c1 * c2 * c3 - s1 * s2 * s3;
					break;

				case 'YXZ':
					this._x = s1 * c2 * c3 + c1 * s2 * s3;
					this._y = c1 * s2 * c3 - s1 * c2 * s3;
					this._z = c1 * c2 * s3 - s1 * s2 * c3;
					this._w = c1 * c2 * c3 + s1 * s2 * s3;
					break;

				case 'ZXY':
					this._x = s1 * c2 * c3 - c1 * s2 * s3;
					this._y = c1 * s2 * c3 + s1 * c2 * s3;
					this._z = c1 * c2 * s3 + s1 * s2 * c3;
					this._w = c1 * c2 * c3 - s1 * s2 * s3;
					break;

				case 'ZYX':
					this._x = s1 * c2 * c3 - c1 * s2 * s3;
					this._y = c1 * s2 * c3 + s1 * c2 * s3;
					this._z = c1 * c2 * s3 - s1 * s2 * c3;
					this._w = c1 * c2 * c3 + s1 * s2 * s3;
					break;

				case 'YZX':
					this._x = s1 * c2 * c3 + c1 * s2 * s3;
					this._y = c1 * s2 * c3 + s1 * c2 * s3;
					this._z = c1 * c2 * s3 - s1 * s2 * c3;
					this._w = c1 * c2 * c3 - s1 * s2 * s3;
					break;

				case 'XZY':
					this._x = s1 * c2 * c3 - c1 * s2 * s3;
					this._y = c1 * s2 * c3 - s1 * c2 * s3;
					this._z = c1 * c2 * s3 + s1 * s2 * c3;
					this._w = c1 * c2 * c3 + s1 * s2 * s3;
					break;

				default:
					console.warn('THREE.Quaternion: .setFromEuler() encountered an unknown order: ' + order);
			}

			if (update !== false) this._onChangeCallback();
			return this;
		};

		_proto.setFromAxisAngle = function setFromAxisAngle(axis, angle) {
			// http://www.euclideanspace.com/maths/geometry/rotations/conversions/angleToQuaternion/index.htm
			// assumes axis is normalized
			var halfAngle = angle / 2,
					s = Math.sin(halfAngle);
			this._x = axis.x * s;
			this._y = axis.y * s;
			this._z = axis.z * s;
			this._w = Math.cos(halfAngle);

			this._onChangeCallback();

			return this;
		};

		_proto.setFromRotationMatrix = function setFromRotationMatrix(m) {
			// http://www.euclideanspace.com/maths/geometry/rotations/conversions/matrixToQuaternion/index.htm
			// assumes the upper 3x3 of m is a pure rotation matrix (i.e, unscaled)
			var te = m.elements,
					m11 = te[0],
					m12 = te[4],
					m13 = te[8],
					m21 = te[1],
					m22 = te[5],
					m23 = te[9],
					m31 = te[2],
					m32 = te[6],
					m33 = te[10],
					trace = m11 + m22 + m33;

			if (trace > 0) {
				var s = 0.5 / Math.sqrt(trace + 1.0);
				this._w = 0.25 / s;
				this._x = (m32 - m23) * s;
				this._y = (m13 - m31) * s;
				this._z = (m21 - m12) * s;
			} else if (m11 > m22 && m11 > m33) {
				var _s = 2.0 * Math.sqrt(1.0 + m11 - m22 - m33);

				this._w = (m32 - m23) / _s;
				this._x = 0.25 * _s;
				this._y = (m12 + m21) / _s;
				this._z = (m13 + m31) / _s;
			} else if (m22 > m33) {
				var _s2 = 2.0 * Math.sqrt(1.0 + m22 - m11 - m33);

				this._w = (m13 - m31) / _s2;
				this._x = (m12 + m21) / _s2;
				this._y = 0.25 * _s2;
				this._z = (m23 + m32) / _s2;
			} else {
				var _s3 = 2.0 * Math.sqrt(1.0 + m33 - m11 - m22);

				this._w = (m21 - m12) / _s3;
				this._x = (m13 + m31) / _s3;
				this._y = (m23 + m32) / _s3;
				this._z = 0.25 * _s3;
			}

			this._onChangeCallback();

			return this;
		};

		_proto.setFromUnitVectors = function setFromUnitVectors(vFrom, vTo) {
			// assumes direction vectors vFrom and vTo are normalized
			var EPS = 0.000001;
			var r = vFrom.dot(vTo) + 1;

			if (r < EPS) {
				r = 0;

				if (Math.abs(vFrom.x) > Math.abs(vFrom.z)) {
					this._x = -vFrom.y;
					this._y = vFrom.x;
					this._z = 0;
					this._w = r;
				} else {
					this._x = 0;
					this._y = -vFrom.z;
					this._z = vFrom.y;
					this._w = r;
				}
			} else {
				// crossVectors( vFrom, vTo ); // inlined to avoid cyclic dependency on Vector3
				this._x = vFrom.y * vTo.z - vFrom.z * vTo.y;
				this._y = vFrom.z * vTo.x - vFrom.x * vTo.z;
				this._z = vFrom.x * vTo.y - vFrom.y * vTo.x;
				this._w = r;
			}

			return this.normalize();
		};

		_proto.angleTo = function angleTo(q) {
			return 2 * Math.acos(Math.abs(MathUtils.clamp(this.dot(q), -1, 1)));
		};

		_proto.rotateTowards = function rotateTowards(q, step) {
			var angle = this.angleTo(q);
			if (angle === 0) return this;
			var t = Math.min(1, step / angle);
			this.slerp(q, t);
			return this;
		};

		_proto.identity = function identity() {
			return this.set(0, 0, 0, 1);
		};

		_proto.inverse = function inverse() {
			// quaternion is assumed to have unit length
			return this.conjugate();
		};

		_proto.conjugate = function conjugate() {
			this._x *= -1;
			this._y *= -1;
			this._z *= -1;

			this._onChangeCallback();

			return this;
		};

		_proto.dot = function dot(v) {
			return this._x * v._x + this._y * v._y + this._z * v._z + this._w * v._w;
		};

		_proto.lengthSq = function lengthSq() {
			return this._x * this._x + this._y * this._y + this._z * this._z + this._w * this._w;
		};

		_proto.length = function length() {
			return Math.sqrt(this._x * this._x + this._y * this._y + this._z * this._z + this._w * this._w);
		};

		_proto.normalize = function normalize() {
			var l = this.length();

			if (l === 0) {
				this._x = 0;
				this._y = 0;
				this._z = 0;
				this._w = 1;
			} else {
				l = 1 / l;
				this._x = this._x * l;
				this._y = this._y * l;
				this._z = this._z * l;
				this._w = this._w * l;
			}

			this._onChangeCallback();

			return this;
		};

		_proto.multiply = function multiply(q, p) {
			if (p !== undefined) {
				console.warn('THREE.Quaternion: .multiply() now only accepts one argument. Use .multiplyQuaternions( a, b ) instead.');
				return this.multiplyQuaternions(q, p);
			}

			return this.multiplyQuaternions(this, q);
		};

		_proto.premultiply = function premultiply(q) {
			return this.multiplyQuaternions(q, this);
		};

		_proto.multiplyQuaternions = function multiplyQuaternions(a, b) {
			// from http://www.euclideanspace.com/maths/algebra/realNormedAlgebra/quaternions/code/index.htm
			var qax = a._x,
					qay = a._y,
					qaz = a._z,
					qaw = a._w;
			var qbx = b._x,
					qby = b._y,
					qbz = b._z,
					qbw = b._w;
			this._x = qax * qbw + qaw * qbx + qay * qbz - qaz * qby;
			this._y = qay * qbw + qaw * qby + qaz * qbx - qax * qbz;
			this._z = qaz * qbw + qaw * qbz + qax * qby - qay * qbx;
			this._w = qaw * qbw - qax * qbx - qay * qby - qaz * qbz;

			this._onChangeCallback();

			return this;
		};

		_proto.slerp = function slerp(qb, t) {
			if (t === 0) return this;
			if (t === 1) return this.copy(qb);
			var x = this._x,
					y = this._y,
					z = this._z,
					w = this._w; // http://www.euclideanspace.com/maths/algebra/realNormedAlgebra/quaternions/slerp/

			var cosHalfTheta = w * qb._w + x * qb._x + y * qb._y + z * qb._z;

			if (cosHalfTheta < 0) {
				this._w = -qb._w;
				this._x = -qb._x;
				this._y = -qb._y;
				this._z = -qb._z;
				cosHalfTheta = -cosHalfTheta;
			} else {
				this.copy(qb);
			}

			if (cosHalfTheta >= 1.0) {
				this._w = w;
				this._x = x;
				this._y = y;
				this._z = z;
				return this;
			}

			var sqrSinHalfTheta = 1.0 - cosHalfTheta * cosHalfTheta;

			if (sqrSinHalfTheta <= Number.EPSILON) {
				var s = 1 - t;
				this._w = s * w + t * this._w;
				this._x = s * x + t * this._x;
				this._y = s * y + t * this._y;
				this._z = s * z + t * this._z;
				this.normalize();

				this._onChangeCallback();

				return this;
			}

			var sinHalfTheta = Math.sqrt(sqrSinHalfTheta);
			var halfTheta = Math.atan2(sinHalfTheta, cosHalfTheta);
			var ratioA = Math.sin((1 - t) * halfTheta) / sinHalfTheta,
					ratioB = Math.sin(t * halfTheta) / sinHalfTheta;
			this._w = w * ratioA + this._w * ratioB;
			this._x = x * ratioA + this._x * ratioB;
			this._y = y * ratioA + this._y * ratioB;
			this._z = z * ratioA + this._z * ratioB;

			this._onChangeCallback();

			return this;
		};

		_proto.equals = function equals(quaternion) {
			return quaternion._x === this._x && quaternion._y === this._y && quaternion._z === this._z && quaternion._w === this._w;
		};

		_proto.fromArray = function fromArray(array, offset) {
			if (offset === undefined) offset = 0;
			this._x = array[offset];
			this._y = array[offset + 1];
			this._z = array[offset + 2];
			this._w = array[offset + 3];

			this._onChangeCallback();

			return this;
		};

		_proto.toArray = function toArray(array, offset) {
			if (array === undefined) array = [];
			if (offset === undefined) offset = 0;
			array[offset] = this._x;
			array[offset + 1] = this._y;
			array[offset + 2] = this._z;
			array[offset + 3] = this._w;
			return array;
		};

		_proto.fromBufferAttribute = function fromBufferAttribute(attribute, index) {
			this._x = attribute.getX(index);
			this._y = attribute.getY(index);
			this._z = attribute.getZ(index);
			this._w = attribute.getW(index);
			return this;
		};

		_proto._onChange = function _onChange(callback) {
			this._onChangeCallback = callback;
			return this;
		};

		_proto._onChangeCallback = function _onChangeCallback() {};

		_createClass(Quaternion, [{
			key: "x",
			get: function get() {
				return this._x;
			},
			set: function set(value) {
				this._x = value;

				this._onChangeCallback();
			}
		}, {
			key: "y",
			get: function get() {
				return this._y;
			},
			set: function set(value) {
				this._y = value;

				this._onChangeCallback();
			}
		}, {
			key: "z",
			get: function get() {
				return this._z;
			},
			set: function set(value) {
				this._z = value;

				this._onChangeCallback();
			}
		}, {
			key: "w",
			get: function get() {
				return this._w;
			},
			set: function set(value) {
				this._w = value;

				this._onChangeCallback();
			}
		}]);

		return Quaternion;
	}();

	var Vector3 = /*#__PURE__*/function () {
		function Vector3(x, y, z) {
			if (x === void 0) {
				x = 0;
			}

			if (y === void 0) {
				y = 0;
			}

			if (z === void 0) {
				z = 0;
			}

			Object.defineProperty(this, 'isVector3', {
				value: true
			});
			this.x = x;
			this.y = y;
			this.z = z;
		}

		var _proto = Vector3.prototype;

		_proto.set = function set(x, y, z) {
			if (z === undefined) z = this.z; // sprite.scale.set(x,y)

			this.x = x;
			this.y = y;
			this.z = z;
			return this;
		};

		_proto.setScalar = function setScalar(scalar) {
			this.x = scalar;
			this.y = scalar;
			this.z = scalar;
			return this;
		};

		_proto.setX = function setX(x) {
			this.x = x;
			return this;
		};

		_proto.setY = function setY(y) {
			this.y = y;
			return this;
		};

		_proto.setZ = function setZ(z) {
			this.z = z;
			return this;
		};

		_proto.setComponent = function setComponent(index, value) {
			switch (index) {
				case 0:
					this.x = value;
					break;

				case 1:
					this.y = value;
					break;

				case 2:
					this.z = value;
					break;

				default:
					throw new Error('index is out of range: ' + index);
			}

			return this;
		};

		_proto.getComponent = function getComponent(index) {
			switch (index) {
				case 0:
					return this.x;

				case 1:
					return this.y;

				case 2:
					return this.z;

				default:
					throw new Error('index is out of range: ' + index);
			}
		};

		_proto.clone = function clone() {
			return new this.constructor(this.x, this.y, this.z);
		};

		_proto.copy = function copy(v) {
			this.x = v.x;
			this.y = v.y;
			this.z = v.z;
			return this;
		};

		_proto.add = function add(v, w) {
			if (w !== undefined) {
				console.warn('THREE.Vector3: .add() now only accepts one argument. Use .addVectors( a, b ) instead.');
				return this.addVectors(v, w);
			}

			this.x += v.x;
			this.y += v.y;
			this.z += v.z;
			return this;
		};

		_proto.addScalar = function addScalar(s) {
			this.x += s;
			this.y += s;
			this.z += s;
			return this;
		};

		_proto.addVectors = function addVectors(a, b) {
			this.x = a.x + b.x;
			this.y = a.y + b.y;
			this.z = a.z + b.z;
			return this;
		};

		_proto.addScaledVector = function addScaledVector(v, s) {
			this.x += v.x * s;
			this.y += v.y * s;
			this.z += v.z * s;
			return this;
		};

		_proto.sub = function sub(v, w) {
			if (w !== undefined) {
				console.warn('THREE.Vector3: .sub() now only accepts one argument. Use .subVectors( a, b ) instead.');
				return this.subVectors(v, w);
			}

			this.x -= v.x;
			this.y -= v.y;
			this.z -= v.z;
			return this;
		};

		_proto.subScalar = function subScalar(s) {
			this.x -= s;
			this.y -= s;
			this.z -= s;
			return this;
		};

		_proto.subVectors = function subVectors(a, b) {
			this.x = a.x - b.x;
			this.y = a.y - b.y;
			this.z = a.z - b.z;
			return this;
		};

		_proto.multiply = function multiply(v, w) {
			if (w !== undefined) {
				console.warn('THREE.Vector3: .multiply() now only accepts one argument. Use .multiplyVectors( a, b ) instead.');
				return this.multiplyVectors(v, w);
			}

			this.x *= v.x;
			this.y *= v.y;
			this.z *= v.z;
			return this;
		};

		_proto.multiplyScalar = function multiplyScalar(scalar) {
			this.x *= scalar;
			this.y *= scalar;
			this.z *= scalar;
			return this;
		};

		_proto.multiplyVectors = function multiplyVectors(a, b) {
			this.x = a.x * b.x;
			this.y = a.y * b.y;
			this.z = a.z * b.z;
			return this;
		};

		_proto.applyEuler = function applyEuler(euler) {
			if (!(euler && euler.isEuler)) {
				console.error('THREE.Vector3: .applyEuler() now expects an Euler rotation rather than a Vector3 and order.');
			}

			return this.applyQuaternion(_quaternion.setFromEuler(euler));
		};

		_proto.applyAxisAngle = function applyAxisAngle(axis, angle) {
			return this.applyQuaternion(_quaternion.setFromAxisAngle(axis, angle));
		};

		_proto.applyMatrix3 = function applyMatrix3(m) {
			var x = this.x,
					y = this.y,
					z = this.z;
			var e = m.elements;
			this.x = e[0] * x + e[3] * y + e[6] * z;
			this.y = e[1] * x + e[4] * y + e[7] * z;
			this.z = e[2] * x + e[5] * y + e[8] * z;
			return this;
		};

		_proto.applyNormalMatrix = function applyNormalMatrix(m) {
			return this.applyMatrix3(m).normalize();
		};

		_proto.applyMatrix4 = function applyMatrix4(m) {
			var x = this.x,
					y = this.y,
					z = this.z;
			var e = m.elements;
			var w = 1 / (e[3] * x + e[7] * y + e[11] * z + e[15]);
			this.x = (e[0] * x + e[4] * y + e[8] * z + e[12]) * w;
			this.y = (e[1] * x + e[5] * y + e[9] * z + e[13]) * w;
			this.z = (e[2] * x + e[6] * y + e[10] * z + e[14]) * w;
			return this;
		};

		_proto.applyQuaternion = function applyQuaternion(q) {
			var x = this.x,
					y = this.y,
					z = this.z;
			var qx = q.x,
					qy = q.y,
					qz = q.z,
					qw = q.w; // calculate quat * vector

			var ix = qw * x + qy * z - qz * y;
			var iy = qw * y + qz * x - qx * z;
			var iz = qw * z + qx * y - qy * x;
			var iw = -qx * x - qy * y - qz * z; // calculate result * inverse quat

			this.x = ix * qw + iw * -qx + iy * -qz - iz * -qy;
			this.y = iy * qw + iw * -qy + iz * -qx - ix * -qz;
			this.z = iz * qw + iw * -qz + ix * -qy - iy * -qx;
			return this;
		};

		_proto.project = function project(camera) {
			return this.applyMatrix4(camera.matrixWorldInverse).applyMatrix4(camera.projectionMatrix);
		};

		_proto.unproject = function unproject(camera) {
			return this.applyMatrix4(camera.projectionMatrixInverse).applyMatrix4(camera.matrixWorld);
		};

		_proto.transformDirection = function transformDirection(m) {
			// input: THREE.Matrix4 affine matrix
			// vector interpreted as a direction
			var x = this.x,
					y = this.y,
					z = this.z;
			var e = m.elements;
			this.x = e[0] * x + e[4] * y + e[8] * z;
			this.y = e[1] * x + e[5] * y + e[9] * z;
			this.z = e[2] * x + e[6] * y + e[10] * z;
			return this.normalize();
		};

		_proto.divide = function divide(v) {
			this.x /= v.x;
			this.y /= v.y;
			this.z /= v.z;
			return this;
		};

		_proto.divideScalar = function divideScalar(scalar) {
			return this.multiplyScalar(1 / scalar);
		};

		_proto.min = function min(v) {
			this.x = Math.min(this.x, v.x);
			this.y = Math.min(this.y, v.y);
			this.z = Math.min(this.z, v.z);
			return this;
		};

		_proto.max = function max(v) {
			this.x = Math.max(this.x, v.x);
			this.y = Math.max(this.y, v.y);
			this.z = Math.max(this.z, v.z);
			return this;
		};

		_proto.clamp = function clamp(min, max) {
			// assumes min < max, componentwise
			this.x = Math.max(min.x, Math.min(max.x, this.x));
			this.y = Math.max(min.y, Math.min(max.y, this.y));
			this.z = Math.max(min.z, Math.min(max.z, this.z));
			return this;
		};

		_proto.clampScalar = function clampScalar(minVal, maxVal) {
			this.x = Math.max(minVal, Math.min(maxVal, this.x));
			this.y = Math.max(minVal, Math.min(maxVal, this.y));
			this.z = Math.max(minVal, Math.min(maxVal, this.z));
			return this;
		};

		_proto.clampLength = function clampLength(min, max) {
			var length = this.length();
			return this.divideScalar(length || 1).multiplyScalar(Math.max(min, Math.min(max, length)));
		};

		_proto.floor = function floor() {
			this.x = Math.floor(this.x);
			this.y = Math.floor(this.y);
			this.z = Math.floor(this.z);
			return this;
		};

		_proto.ceil = function ceil() {
			this.x = Math.ceil(this.x);
			this.y = Math.ceil(this.y);
			this.z = Math.ceil(this.z);
			return this;
		};

		_proto.round = function round() {
			this.x = Math.round(this.x);
			this.y = Math.round(this.y);
			this.z = Math.round(this.z);
			return this;
		};

		_proto.roundToZero = function roundToZero() {
			this.x = this.x < 0 ? Math.ceil(this.x) : Math.floor(this.x);
			this.y = this.y < 0 ? Math.ceil(this.y) : Math.floor(this.y);
			this.z = this.z < 0 ? Math.ceil(this.z) : Math.floor(this.z);
			return this;
		};

		_proto.negate = function negate() {
			this.x = -this.x;
			this.y = -this.y;
			this.z = -this.z;
			return this;
		};

		_proto.dot = function dot(v) {
			return this.x * v.x + this.y * v.y + this.z * v.z;
		} // TODO lengthSquared?
		;

		_proto.lengthSq = function lengthSq() {
			return this.x * this.x + this.y * this.y + this.z * this.z;
		};

		_proto.length = function length() {
			return Math.sqrt(this.x * this.x + this.y * this.y + this.z * this.z);
		};

		_proto.manhattanLength = function manhattanLength() {
			return Math.abs(this.x) + Math.abs(this.y) + Math.abs(this.z);
		};

		_proto.normalize = function normalize() {
			return this.divideScalar(this.length() || 1);
		};

		_proto.setLength = function setLength(length) {
			return this.normalize().multiplyScalar(length);
		};

		_proto.lerp = function lerp(v, alpha) {
			this.x += (v.x - this.x) * alpha;
			this.y += (v.y - this.y) * alpha;
			this.z += (v.z - this.z) * alpha;
			return this;
		};

		_proto.lerpVectors = function lerpVectors(v1, v2, alpha) {
			this.x = v1.x + (v2.x - v1.x) * alpha;
			this.y = v1.y + (v2.y - v1.y) * alpha;
			this.z = v1.z + (v2.z - v1.z) * alpha;
			return this;
		};

		_proto.cross = function cross(v, w) {
			if (w !== undefined) {
				console.warn('THREE.Vector3: .cross() now only accepts one argument. Use .crossVectors( a, b ) instead.');
				return this.crossVectors(v, w);
			}

			return this.crossVectors(this, v);
		};

		_proto.crossVectors = function crossVectors(a, b) {
			var ax = a.x,
					ay = a.y,
					az = a.z;
			var bx = b.x,
					by = b.y,
					bz = b.z;
			this.x = ay * bz - az * by;
			this.y = az * bx - ax * bz;
			this.z = ax * by - ay * bx;
			return this;
		};

		_proto.projectOnVector = function projectOnVector(v) {
			var denominator = v.lengthSq();
			if (denominator === 0) return this.set(0, 0, 0);
			var scalar = v.dot(this) / denominator;
			return this.copy(v).multiplyScalar(scalar);
		};

		_proto.projectOnPlane = function projectOnPlane(planeNormal) {
			_vector.copy(this).projectOnVector(planeNormal);

			return this.sub(_vector);
		};

		_proto.reflect = function reflect(normal) {
			// reflect incident vector off plane orthogonal to normal
			// normal is assumed to have unit length
			return this.sub(_vector.copy(normal).multiplyScalar(2 * this.dot(normal)));
		};

		_proto.angleTo = function angleTo(v) {
			var denominator = Math.sqrt(this.lengthSq() * v.lengthSq());
			if (denominator === 0) return Math.PI / 2;
			var theta = this.dot(v) / denominator; // clamp, to handle numerical problems

			return Math.acos(MathUtils.clamp(theta, -1, 1));
		};

		_proto.distanceTo = function distanceTo(v) {
			return Math.sqrt(this.distanceToSquared(v));
		};

		_proto.distanceToSquared = function distanceToSquared(v) {
			var dx = this.x - v.x,
					dy = this.y - v.y,
					dz = this.z - v.z;
			return dx * dx + dy * dy + dz * dz;
		};

		_proto.manhattanDistanceTo = function manhattanDistanceTo(v) {
			return Math.abs(this.x - v.x) + Math.abs(this.y - v.y) + Math.abs(this.z - v.z);
		};

		_proto.setFromSpherical = function setFromSpherical(s) {
			return this.setFromSphericalCoords(s.radius, s.phi, s.theta);
		};

		_proto.setFromSphericalCoords = function setFromSphericalCoords(radius, phi, theta) {
			var sinPhiRadius = Math.sin(phi) * radius;
			this.x = sinPhiRadius * Math.sin(theta);
			this.y = Math.cos(phi) * radius;
			this.z = sinPhiRadius * Math.cos(theta);
			return this;
		};

		_proto.setFromCylindrical = function setFromCylindrical(c) {
			return this.setFromCylindricalCoords(c.radius, c.theta, c.y);
		};

		_proto.setFromCylindricalCoords = function setFromCylindricalCoords(radius, theta, y) {
			this.x = radius * Math.sin(theta);
			this.y = y;
			this.z = radius * Math.cos(theta);
			return this;
		};

		_proto.setFromMatrixPosition = function setFromMatrixPosition(m) {
			var e = m.elements;
			this.x = e[12];
			this.y = e[13];
			this.z = e[14];
			return this;
		};

		_proto.setFromMatrixScale = function setFromMatrixScale(m) {
			var sx = this.setFromMatrixColumn(m, 0).length();
			var sy = this.setFromMatrixColumn(m, 1).length();
			var sz = this.setFromMatrixColumn(m, 2).length();
			this.x = sx;
			this.y = sy;
			this.z = sz;
			return this;
		};

		_proto.setFromMatrixColumn = function setFromMatrixColumn(m, index) {
			return this.fromArray(m.elements, index * 4);
		};

		_proto.setFromMatrix3Column = function setFromMatrix3Column(m, index) {
			return this.fromArray(m.elements, index * 3);
		};

		_proto.equals = function equals(v) {
			return v.x === this.x && v.y === this.y && v.z === this.z;
		};

		_proto.fromArray = function fromArray(array, offset) {
			if (offset === undefined) offset = 0;
			this.x = array[offset];
			this.y = array[offset + 1];
			this.z = array[offset + 2];
			return this;
		};

		_proto.toArray = function toArray(array, offset) {
			if (array === undefined) array = [];
			if (offset === undefined) offset = 0;
			array[offset] = this.x;
			array[offset + 1] = this.y;
			array[offset + 2] = this.z;
			return array;
		};

		_proto.fromBufferAttribute = function fromBufferAttribute(attribute, index, offset) {
			if (offset !== undefined) {
				console.warn('THREE.Vector3: offset has been removed from .fromBufferAttribute().');
			}

			this.x = attribute.getX(index);
			this.y = attribute.getY(index);
			this.z = attribute.getZ(index);
			return this;
		};

		_proto.random = function random() {
			this.x = Math.random();
			this.y = Math.random();
			this.z = Math.random();
			return this;
		};

		return Vector3;
	}();

	var _vector = /*@__PURE__*/new Vector3();

	var _quaternion = /*@__PURE__*/new Quaternion();

	var Matrix4 = /*#__PURE__*/function () {
		function Matrix4() {
			Object.defineProperty(this, 'isMatrix4', {
				value: true
			});
			this.elements = [1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1];

			if (arguments.length > 0) {
				console.error('THREE.Matrix4: the constructor no longer reads arguments. use .set() instead.');
			}
		}

		var _proto = Matrix4.prototype;

		_proto.set = function set(n11, n12, n13, n14, n21, n22, n23, n24, n31, n32, n33, n34, n41, n42, n43, n44) {
			var te = this.elements;
			te[0] = n11;
			te[4] = n12;
			te[8] = n13;
			te[12] = n14;
			te[1] = n21;
			te[5] = n22;
			te[9] = n23;
			te[13] = n24;
			te[2] = n31;
			te[6] = n32;
			te[10] = n33;
			te[14] = n34;
			te[3] = n41;
			te[7] = n42;
			te[11] = n43;
			te[15] = n44;
			return this;
		};

		_proto.identity = function identity() {
			this.set(1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1);
			return this;
		};

		_proto.clone = function clone() {
			return new Matrix4().fromArray(this.elements);
		};

		_proto.copy = function copy(m) {
			var te = this.elements;
			var me = m.elements;
			te[0] = me[0];
			te[1] = me[1];
			te[2] = me[2];
			te[3] = me[3];
			te[4] = me[4];
			te[5] = me[5];
			te[6] = me[6];
			te[7] = me[7];
			te[8] = me[8];
			te[9] = me[9];
			te[10] = me[10];
			te[11] = me[11];
			te[12] = me[12];
			te[13] = me[13];
			te[14] = me[14];
			te[15] = me[15];
			return this;
		};

		_proto.copyPosition = function copyPosition(m) {
			var te = this.elements,
					me = m.elements;
			te[12] = me[12];
			te[13] = me[13];
			te[14] = me[14];
			return this;
		};

		_proto.extractBasis = function extractBasis(xAxis, yAxis, zAxis) {
			xAxis.setFromMatrixColumn(this, 0);
			yAxis.setFromMatrixColumn(this, 1);
			zAxis.setFromMatrixColumn(this, 2);
			return this;
		};

		_proto.makeBasis = function makeBasis(xAxis, yAxis, zAxis) {
			this.set(xAxis.x, yAxis.x, zAxis.x, 0, xAxis.y, yAxis.y, zAxis.y, 0, xAxis.z, yAxis.z, zAxis.z, 0, 0, 0, 0, 1);
			return this;
		};

		_proto.extractRotation = function extractRotation(m) {
			// this method does not support reflection matrices
			var te = this.elements;
			var me = m.elements;

			var scaleX = 1 / _v1.setFromMatrixColumn(m, 0).length();

			var scaleY = 1 / _v1.setFromMatrixColumn(m, 1).length();

			var scaleZ = 1 / _v1.setFromMatrixColumn(m, 2).length();

			te[0] = me[0] * scaleX;
			te[1] = me[1] * scaleX;
			te[2] = me[2] * scaleX;
			te[3] = 0;
			te[4] = me[4] * scaleY;
			te[5] = me[5] * scaleY;
			te[6] = me[6] * scaleY;
			te[7] = 0;
			te[8] = me[8] * scaleZ;
			te[9] = me[9] * scaleZ;
			te[10] = me[10] * scaleZ;
			te[11] = 0;
			te[12] = 0;
			te[13] = 0;
			te[14] = 0;
			te[15] = 1;
			return this;
		};

		_proto.makeRotationFromEuler = function makeRotationFromEuler(euler) {
			if (!(euler && euler.isEuler)) {
				console.error('THREE.Matrix4: .makeRotationFromEuler() now expects a Euler rotation rather than a Vector3 and order.');
			}

			var te = this.elements;
			var x = euler.x,
					y = euler.y,
					z = euler.z;
			var a = Math.cos(x),
					b = Math.sin(x);
			var c = Math.cos(y),
					d = Math.sin(y);
			var e = Math.cos(z),
					f = Math.sin(z);

			if (euler.order === 'XYZ') {
				var ae = a * e,
						af = a * f,
						be = b * e,
						bf = b * f;
				te[0] = c * e;
				te[4] = -c * f;
				te[8] = d;
				te[1] = af + be * d;
				te[5] = ae - bf * d;
				te[9] = -b * c;
				te[2] = bf - ae * d;
				te[6] = be + af * d;
				te[10] = a * c;
			} else if (euler.order === 'YXZ') {
				var ce = c * e,
						cf = c * f,
						de = d * e,
						df = d * f;
				te[0] = ce + df * b;
				te[4] = de * b - cf;
				te[8] = a * d;
				te[1] = a * f;
				te[5] = a * e;
				te[9] = -b;
				te[2] = cf * b - de;
				te[6] = df + ce * b;
				te[10] = a * c;
			} else if (euler.order === 'ZXY') {
				var _ce = c * e,
						_cf = c * f,
						_de = d * e,
						_df = d * f;

				te[0] = _ce - _df * b;
				te[4] = -a * f;
				te[8] = _de + _cf * b;
				te[1] = _cf + _de * b;
				te[5] = a * e;
				te[9] = _df - _ce * b;
				te[2] = -a * d;
				te[6] = b;
				te[10] = a * c;
			} else if (euler.order === 'ZYX') {
				var _ae = a * e,
						_af = a * f,
						_be = b * e,
						_bf = b * f;

				te[0] = c * e;
				te[4] = _be * d - _af;
				te[8] = _ae * d + _bf;
				te[1] = c * f;
				te[5] = _bf * d + _ae;
				te[9] = _af * d - _be;
				te[2] = -d;
				te[6] = b * c;
				te[10] = a * c;
			} else if (euler.order === 'YZX') {
				var ac = a * c,
						ad = a * d,
						bc = b * c,
						bd = b * d;
				te[0] = c * e;
				te[4] = bd - ac * f;
				te[8] = bc * f + ad;
				te[1] = f;
				te[5] = a * e;
				te[9] = -b * e;
				te[2] = -d * e;
				te[6] = ad * f + bc;
				te[10] = ac - bd * f;
			} else if (euler.order === 'XZY') {
				var _ac = a * c,
						_ad = a * d,
						_bc = b * c,
						_bd = b * d;

				te[0] = c * e;
				te[4] = -f;
				te[8] = d * e;
				te[1] = _ac * f + _bd;
				te[5] = a * e;
				te[9] = _ad * f - _bc;
				te[2] = _bc * f - _ad;
				te[6] = b * e;
				te[10] = _bd * f + _ac;
			} // bottom row


			te[3] = 0;
			te[7] = 0;
			te[11] = 0; // last column

			te[12] = 0;
			te[13] = 0;
			te[14] = 0;
			te[15] = 1;
			return this;
		};

		_proto.makeRotationFromQuaternion = function makeRotationFromQuaternion(q) {
			return this.compose(_zero, q, _one);
		};

		_proto.lookAt = function lookAt(eye, target, up) {
			var te = this.elements;

			_z.subVectors(eye, target);

			if (_z.lengthSq() === 0) {
				// eye and target are in the same position
				_z.z = 1;
			}

			_z.normalize();

			_x.crossVectors(up, _z);

			if (_x.lengthSq() === 0) {
				// up and z are parallel
				if (Math.abs(up.z) === 1) {
					_z.x += 0.0001;
				} else {
					_z.z += 0.0001;
				}

				_z.normalize();

				_x.crossVectors(up, _z);
			}

			_x.normalize();

			_y.crossVectors(_z, _x);

			te[0] = _x.x;
			te[4] = _y.x;
			te[8] = _z.x;
			te[1] = _x.y;
			te[5] = _y.y;
			te[9] = _z.y;
			te[2] = _x.z;
			te[6] = _y.z;
			te[10] = _z.z;
			return this;
		};

		_proto.multiply = function multiply(m, n) {
			if (n !== undefined) {
				console.warn('THREE.Matrix4: .multiply() now only accepts one argument. Use .multiplyMatrices( a, b ) instead.');
				return this.multiplyMatrices(m, n);
			}

			return this.multiplyMatrices(this, m);
		};

		_proto.premultiply = function premultiply(m) {
			return this.multiplyMatrices(m, this);
		};

		_proto.multiplyMatrices = function multiplyMatrices(a, b) {
			var ae = a.elements;
			var be = b.elements;
			var te = this.elements;
			var a11 = ae[0],
					a12 = ae[4],
					a13 = ae[8],
					a14 = ae[12];
			var a21 = ae[1],
					a22 = ae[5],
					a23 = ae[9],
					a24 = ae[13];
			var a31 = ae[2],
					a32 = ae[6],
					a33 = ae[10],
					a34 = ae[14];
			var a41 = ae[3],
					a42 = ae[7],
					a43 = ae[11],
					a44 = ae[15];
			var b11 = be[0],
					b12 = be[4],
					b13 = be[8],
					b14 = be[12];
			var b21 = be[1],
					b22 = be[5],
					b23 = be[9],
					b24 = be[13];
			var b31 = be[2],
					b32 = be[6],
					b33 = be[10],
					b34 = be[14];
			var b41 = be[3],
					b42 = be[7],
					b43 = be[11],
					b44 = be[15];
			te[0] = a11 * b11 + a12 * b21 + a13 * b31 + a14 * b41;
			te[4] = a11 * b12 + a12 * b22 + a13 * b32 + a14 * b42;
			te[8] = a11 * b13 + a12 * b23 + a13 * b33 + a14 * b43;
			te[12] = a11 * b14 + a12 * b24 + a13 * b34 + a14 * b44;
			te[1] = a21 * b11 + a22 * b21 + a23 * b31 + a24 * b41;
			te[5] = a21 * b12 + a22 * b22 + a23 * b32 + a24 * b42;
			te[9] = a21 * b13 + a22 * b23 + a23 * b33 + a24 * b43;
			te[13] = a21 * b14 + a22 * b24 + a23 * b34 + a24 * b44;
			te[2] = a31 * b11 + a32 * b21 + a33 * b31 + a34 * b41;
			te[6] = a31 * b12 + a32 * b22 + a33 * b32 + a34 * b42;
			te[10] = a31 * b13 + a32 * b23 + a33 * b33 + a34 * b43;
			te[14] = a31 * b14 + a32 * b24 + a33 * b34 + a34 * b44;
			te[3] = a41 * b11 + a42 * b21 + a43 * b31 + a44 * b41;
			te[7] = a41 * b12 + a42 * b22 + a43 * b32 + a44 * b42;
			te[11] = a41 * b13 + a42 * b23 + a43 * b33 + a44 * b43;
			te[15] = a41 * b14 + a42 * b24 + a43 * b34 + a44 * b44;
			return this;
		};

		_proto.multiplyScalar = function multiplyScalar(s) {
			var te = this.elements;
			te[0] *= s;
			te[4] *= s;
			te[8] *= s;
			te[12] *= s;
			te[1] *= s;
			te[5] *= s;
			te[9] *= s;
			te[13] *= s;
			te[2] *= s;
			te[6] *= s;
			te[10] *= s;
			te[14] *= s;
			te[3] *= s;
			te[7] *= s;
			te[11] *= s;
			te[15] *= s;
			return this;
		};

		_proto.determinant = function determinant() {
			var te = this.elements;
			var n11 = te[0],
					n12 = te[4],
					n13 = te[8],
					n14 = te[12];
			var n21 = te[1],
					n22 = te[5],
					n23 = te[9],
					n24 = te[13];
			var n31 = te[2],
					n32 = te[6],
					n33 = te[10],
					n34 = te[14];
			var n41 = te[3],
					n42 = te[7],
					n43 = te[11],
					n44 = te[15]; //TODO: make this more efficient
			//( based on http://www.euclideanspace.com/maths/algebra/matrix/functions/inverse/fourD/index.htm )

			return n41 * (+n14 * n23 * n32 - n13 * n24 * n32 - n14 * n22 * n33 + n12 * n24 * n33 + n13 * n22 * n34 - n12 * n23 * n34) + n42 * (+n11 * n23 * n34 - n11 * n24 * n33 + n14 * n21 * n33 - n13 * n21 * n34 + n13 * n24 * n31 - n14 * n23 * n31) + n43 * (+n11 * n24 * n32 - n11 * n22 * n34 - n14 * n21 * n32 + n12 * n21 * n34 + n14 * n22 * n31 - n12 * n24 * n31) + n44 * (-n13 * n22 * n31 - n11 * n23 * n32 + n11 * n22 * n33 + n13 * n21 * n32 - n12 * n21 * n33 + n12 * n23 * n31);
		};

		_proto.transpose = function transpose() {
			var te = this.elements;
			var tmp;
			tmp = te[1];
			te[1] = te[4];
			te[4] = tmp;
			tmp = te[2];
			te[2] = te[8];
			te[8] = tmp;
			tmp = te[6];
			te[6] = te[9];
			te[9] = tmp;
			tmp = te[3];
			te[3] = te[12];
			te[12] = tmp;
			tmp = te[7];
			te[7] = te[13];
			te[13] = tmp;
			tmp = te[11];
			te[11] = te[14];
			te[14] = tmp;
			return this;
		};

		_proto.setPosition = function setPosition(x, y, z) {
			var te = this.elements;

			if (x.isVector3) {
				te[12] = x.x;
				te[13] = x.y;
				te[14] = x.z;
			} else {
				te[12] = x;
				te[13] = y;
				te[14] = z;
			}

			return this;
		};

		_proto.getInverse = function getInverse(m, throwOnDegenerate) {
			if (throwOnDegenerate !== undefined) {
				console.warn("THREE.Matrix4: .getInverse() can no longer be configured to throw on degenerate.");
			} // based on http://www.euclideanspace.com/maths/algebra/matrix/functions/inverse/fourD/index.htm


			var te = this.elements,
					me = m.elements,
					n11 = me[0],
					n21 = me[1],
					n31 = me[2],
					n41 = me[3],
					n12 = me[4],
					n22 = me[5],
					n32 = me[6],
					n42 = me[7],
					n13 = me[8],
					n23 = me[9],
					n33 = me[10],
					n43 = me[11],
					n14 = me[12],
					n24 = me[13],
					n34 = me[14],
					n44 = me[15],
					t11 = n23 * n34 * n42 - n24 * n33 * n42 + n24 * n32 * n43 - n22 * n34 * n43 - n23 * n32 * n44 + n22 * n33 * n44,
					t12 = n14 * n33 * n42 - n13 * n34 * n42 - n14 * n32 * n43 + n12 * n34 * n43 + n13 * n32 * n44 - n12 * n33 * n44,
					t13 = n13 * n24 * n42 - n14 * n23 * n42 + n14 * n22 * n43 - n12 * n24 * n43 - n13 * n22 * n44 + n12 * n23 * n44,
					t14 = n14 * n23 * n32 - n13 * n24 * n32 - n14 * n22 * n33 + n12 * n24 * n33 + n13 * n22 * n34 - n12 * n23 * n34;
			var det = n11 * t11 + n21 * t12 + n31 * t13 + n41 * t14;
			if (det === 0) return this.set(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);
			var detInv = 1 / det;
			te[0] = t11 * detInv;
			te[1] = (n24 * n33 * n41 - n23 * n34 * n41 - n24 * n31 * n43 + n21 * n34 * n43 + n23 * n31 * n44 - n21 * n33 * n44) * detInv;
			te[2] = (n22 * n34 * n41 - n24 * n32 * n41 + n24 * n31 * n42 - n21 * n34 * n42 - n22 * n31 * n44 + n21 * n32 * n44) * detInv;
			te[3] = (n23 * n32 * n41 - n22 * n33 * n41 - n23 * n31 * n42 + n21 * n33 * n42 + n22 * n31 * n43 - n21 * n32 * n43) * detInv;
			te[4] = t12 * detInv;
			te[5] = (n13 * n34 * n41 - n14 * n33 * n41 + n14 * n31 * n43 - n11 * n34 * n43 - n13 * n31 * n44 + n11 * n33 * n44) * detInv;
			te[6] = (n14 * n32 * n41 - n12 * n34 * n41 - n14 * n31 * n42 + n11 * n34 * n42 + n12 * n31 * n44 - n11 * n32 * n44) * detInv;
			te[7] = (n12 * n33 * n41 - n13 * n32 * n41 + n13 * n31 * n42 - n11 * n33 * n42 - n12 * n31 * n43 + n11 * n32 * n43) * detInv;
			te[8] = t13 * detInv;
			te[9] = (n14 * n23 * n41 - n13 * n24 * n41 - n14 * n21 * n43 + n11 * n24 * n43 + n13 * n21 * n44 - n11 * n23 * n44) * detInv;
			te[10] = (n12 * n24 * n41 - n14 * n22 * n41 + n14 * n21 * n42 - n11 * n24 * n42 - n12 * n21 * n44 + n11 * n22 * n44) * detInv;
			te[11] = (n13 * n22 * n41 - n12 * n23 * n41 - n13 * n21 * n42 + n11 * n23 * n42 + n12 * n21 * n43 - n11 * n22 * n43) * detInv;
			te[12] = t14 * detInv;
			te[13] = (n13 * n24 * n31 - n14 * n23 * n31 + n14 * n21 * n33 - n11 * n24 * n33 - n13 * n21 * n34 + n11 * n23 * n34) * detInv;
			te[14] = (n14 * n22 * n31 - n12 * n24 * n31 - n14 * n21 * n32 + n11 * n24 * n32 + n12 * n21 * n34 - n11 * n22 * n34) * detInv;
			te[15] = (n12 * n23 * n31 - n13 * n22 * n31 + n13 * n21 * n32 - n11 * n23 * n32 - n12 * n21 * n33 + n11 * n22 * n33) * detInv;
			return this;
		};

		_proto.scale = function scale(v) {
			var te = this.elements;
			var x = v.x,
					y = v.y,
					z = v.z;
			te[0] *= x;
			te[4] *= y;
			te[8] *= z;
			te[1] *= x;
			te[5] *= y;
			te[9] *= z;
			te[2] *= x;
			te[6] *= y;
			te[10] *= z;
			te[3] *= x;
			te[7] *= y;
			te[11] *= z;
			return this;
		};

		_proto.getMaxScaleOnAxis = function getMaxScaleOnAxis() {
			var te = this.elements;
			var scaleXSq = te[0] * te[0] + te[1] * te[1] + te[2] * te[2];
			var scaleYSq = te[4] * te[4] + te[5] * te[5] + te[6] * te[6];
			var scaleZSq = te[8] * te[8] + te[9] * te[9] + te[10] * te[10];
			return Math.sqrt(Math.max(scaleXSq, scaleYSq, scaleZSq));
		};

		_proto.makeTranslation = function makeTranslation(x, y, z) {
			this.set(1, 0, 0, x, 0, 1, 0, y, 0, 0, 1, z, 0, 0, 0, 1);
			return this;
		};

		_proto.makeRotationX = function makeRotationX(theta) {
			var c = Math.cos(theta),
					s = Math.sin(theta);
			this.set(1, 0, 0, 0, 0, c, -s, 0, 0, s, c, 0, 0, 0, 0, 1);
			return this;
		};

		_proto.makeRotationY = function makeRotationY(theta) {
			var c = Math.cos(theta),
					s = Math.sin(theta);
			this.set(c, 0, s, 0, 0, 1, 0, 0, -s, 0, c, 0, 0, 0, 0, 1);
			return this;
		};

		_proto.makeRotationZ = function makeRotationZ(theta) {
			var c = Math.cos(theta),
					s = Math.sin(theta);
			this.set(c, -s, 0, 0, s, c, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1);
			return this;
		};

		_proto.makeRotationAxis = function makeRotationAxis(axis, angle) {
			// Based on http://www.gamedev.net/reference/articles/article1199.asp
			var c = Math.cos(angle);
			var s = Math.sin(angle);
			var t = 1 - c;
			var x = axis.x,
					y = axis.y,
					z = axis.z;
			var tx = t * x,
					ty = t * y;
			this.set(tx * x + c, tx * y - s * z, tx * z + s * y, 0, tx * y + s * z, ty * y + c, ty * z - s * x, 0, tx * z - s * y, ty * z + s * x, t * z * z + c, 0, 0, 0, 0, 1);
			return this;
		};

		_proto.makeScale = function makeScale(x, y, z) {
			this.set(x, 0, 0, 0, 0, y, 0, 0, 0, 0, z, 0, 0, 0, 0, 1);
			return this;
		};

		_proto.makeShear = function makeShear(x, y, z) {
			this.set(1, y, z, 0, x, 1, z, 0, x, y, 1, 0, 0, 0, 0, 1);
			return this;
		};

		_proto.compose = function compose(position, quaternion, scale) {
			var te = this.elements;
			var x = quaternion._x,
					y = quaternion._y,
					z = quaternion._z,
					w = quaternion._w;
			var x2 = x + x,
					y2 = y + y,
					z2 = z + z;
			var xx = x * x2,
					xy = x * y2,
					xz = x * z2;
			var yy = y * y2,
					yz = y * z2,
					zz = z * z2;
			var wx = w * x2,
					wy = w * y2,
					wz = w * z2;
			var sx = scale.x,
					sy = scale.y,
					sz = scale.z;
			te[0] = (1 - (yy + zz)) * sx;
			te[1] = (xy + wz) * sx;
			te[2] = (xz - wy) * sx;
			te[3] = 0;
			te[4] = (xy - wz) * sy;
			te[5] = (1 - (xx + zz)) * sy;
			te[6] = (yz + wx) * sy;
			te[7] = 0;
			te[8] = (xz + wy) * sz;
			te[9] = (yz - wx) * sz;
			te[10] = (1 - (xx + yy)) * sz;
			te[11] = 0;
			te[12] = position.x;
			te[13] = position.y;
			te[14] = position.z;
			te[15] = 1;
			return this;
		};

		_proto.decompose = function decompose(position, quaternion, scale) {
			var te = this.elements;

			var sx = _v1.set(te[0], te[1], te[2]).length();

			var sy = _v1.set(te[4], te[5], te[6]).length();

			var sz = _v1.set(te[8], te[9], te[10]).length(); // if determine is negative, we need to invert one scale


			var det = this.determinant();
			if (det < 0) sx = -sx;
			position.x = te[12];
			position.y = te[13];
			position.z = te[14]; // scale the rotation part

			_m1.copy(this);

			var invSX = 1 / sx;
			var invSY = 1 / sy;
			var invSZ = 1 / sz;
			_m1.elements[0] *= invSX;
			_m1.elements[1] *= invSX;
			_m1.elements[2] *= invSX;
			_m1.elements[4] *= invSY;
			_m1.elements[5] *= invSY;
			_m1.elements[6] *= invSY;
			_m1.elements[8] *= invSZ;
			_m1.elements[9] *= invSZ;
			_m1.elements[10] *= invSZ;
			quaternion.setFromRotationMatrix(_m1);
			scale.x = sx;
			scale.y = sy;
			scale.z = sz;
			return this;
		};

		_proto.makePerspective = function makePerspective(left, right, top, bottom, near, far) {
			if (far === undefined) {
				console.warn('THREE.Matrix4: .makePerspective() has been redefined and has a new signature. Please check the docs.');
			}

			var te = this.elements;
			var x = 2 * near / (right - left);
			var y = 2 * near / (top - bottom);
			var a = (right + left) / (right - left);
			var b = (top + bottom) / (top - bottom);
			var c = -(far + near) / (far - near);
			var d = -2 * far * near / (far - near);
			te[0] = x;
			te[4] = 0;
			te[8] = a;
			te[12] = 0;
			te[1] = 0;
			te[5] = y;
			te[9] = b;
			te[13] = 0;
			te[2] = 0;
			te[6] = 0;
			te[10] = c;
			te[14] = d;
			te[3] = 0;
			te[7] = 0;
			te[11] = -1;
			te[15] = 0;
			return this;
		};

		_proto.makeOrthographic = function makeOrthographic(left, right, top, bottom, near, far) {
			var te = this.elements;
			var w = 1.0 / (right - left);
			var h = 1.0 / (top - bottom);
			var p = 1.0 / (far - near);
			var x = (right + left) * w;
			var y = (top + bottom) * h;
			var z = (far + near) * p;
			te[0] = 2 * w;
			te[4] = 0;
			te[8] = 0;
			te[12] = -x;
			te[1] = 0;
			te[5] = 2 * h;
			te[9] = 0;
			te[13] = -y;
			te[2] = 0;
			te[6] = 0;
			te[10] = -2 * p;
			te[14] = -z;
			te[3] = 0;
			te[7] = 0;
			te[11] = 0;
			te[15] = 1;
			return this;
		};

		_proto.equals = function equals(matrix) {
			var te = this.elements;
			var me = matrix.elements;

			for (var i = 0; i < 16; i++) {
				if (te[i] !== me[i]) return false;
			}

			return true;
		};

		_proto.fromArray = function fromArray(array, offset) {
			if (offset === undefined) offset = 0;

			for (var i = 0; i < 16; i++) {
				this.elements[i] = array[i + offset];
			}

			return this;
		};

		_proto.toArray = function toArray(array, offset) {
			if (array === undefined) array = [];
			if (offset === undefined) offset = 0;
			var te = this.elements;
			array[offset] = te[0];
			array[offset + 1] = te[1];
			array[offset + 2] = te[2];
			array[offset + 3] = te[3];
			array[offset + 4] = te[4];
			array[offset + 5] = te[5];
			array[offset + 6] = te[6];
			array[offset + 7] = te[7];
			array[offset + 8] = te[8];
			array[offset + 9] = te[9];
			array[offset + 10] = te[10];
			array[offset + 11] = te[11];
			array[offset + 12] = te[12];
			array[offset + 13] = te[13];
			array[offset + 14] = te[14];
			array[offset + 15] = te[15];
			return array;
		};

		return Matrix4;
	}();

	var _v1 = /*@__PURE__*/new Vector3();

	var _m1 = /*@__PURE__*/new Matrix4();

	var _zero = /*@__PURE__*/new Vector3(0, 0, 0);

	var _one = /*@__PURE__*/new Vector3(1, 1, 1);

	var _x = /*@__PURE__*/new Vector3();

	var _y = /*@__PURE__*/new Vector3();

	var _z = /*@__PURE__*/new Vector3();

	/**
	 * Extensible curve object.
	 *
	 * Some common of curve methods:
	 * .getPoint( t, optionalTarget ), .getTangent( t, optionalTarget )
	 * .getPointAt( u, optionalTarget ), .getTangentAt( u, optionalTarget )
	 * .getPoints(), .getSpacedPoints()
	 * .getLength()
	 * .updateArcLengths()
	 *
	 * This following curves inherit from THREE.Curve:
	 *
	 * -- 2D curves --
	 * THREE.ArcCurve
	 * THREE.CubicBezierCurve
	 * THREE.EllipseCurve
	 * THREE.LineCurve
	 * THREE.QuadraticBezierCurve
	 * THREE.SplineCurve
	 *
	 * -- 3D curves --
	 * THREE.CatmullRomCurve3
	 * THREE.CubicBezierCurve3
	 * THREE.LineCurve3
	 * THREE.QuadraticBezierCurve3
	 *
	 * A series of curves can be represented as a THREE.CurvePath.
	 *
	 **/

	function Curve() {
		this.type = 'Curve';
		this.arcLengthDivisions = 200;
	}

	Object.assign(Curve.prototype, {
		// Virtual base class method to overwrite and implement in subclasses
		//	- t [0 .. 1]
		getPoint: function getPoint()
		/* t, optionalTarget */
		{
			console.warn('THREE.Curve: .getPoint() not implemented.');
			return null;
		},
		// Get point at relative position in curve according to arc length
		// - u [0 .. 1]
		getPointAt: function getPointAt(u, optionalTarget) {
			var t = this.getUtoTmapping(u);
			return this.getPoint(t, optionalTarget);
		},
		// Get sequence of points using getPoint( t )
		getPoints: function getPoints(divisions) {
			if (divisions === undefined) divisions = 5;
			var points = [];

			for (var d = 0; d <= divisions; d++) {
				points.push(this.getPoint(d / divisions));
			}

			return points;
		},
		// Get sequence of points using getPointAt( u )
		getSpacedPoints: function getSpacedPoints(divisions) {
			if (divisions === undefined) divisions = 5;
			var points = [];

			for (var d = 0; d <= divisions; d++) {
				points.push(this.getPointAt(d / divisions));
			}

			return points;
		},
		// Get total curve arc length
		getLength: function getLength() {
			var lengths = this.getLengths();
			return lengths[lengths.length - 1];
		},
		// Get list of cumulative segment lengths
		getLengths: function getLengths(divisions) {
			if (divisions === undefined) divisions = this.arcLengthDivisions;

			if (this.cacheArcLengths && this.cacheArcLengths.length === divisions + 1 && !this.needsUpdate) {
				return this.cacheArcLengths;
			}

			this.needsUpdate = false;
			var cache = [];
			var current,
					last = this.getPoint(0);
			var sum = 0;
			cache.push(0);

			for (var p = 1; p <= divisions; p++) {
				current = this.getPoint(p / divisions);
				sum += current.distanceTo(last);
				cache.push(sum);
				last = current;
			}

			this.cacheArcLengths = cache;
			return cache; // { sums: cache, sum: sum }; Sum is in the last element.
		},
		updateArcLengths: function updateArcLengths() {
			this.needsUpdate = true;
			this.getLengths();
		},
		// Given u ( 0 .. 1 ), get a t to find p. This gives you points which are equidistant
		getUtoTmapping: function getUtoTmapping(u, distance) {
			var arcLengths = this.getLengths();
			var i = 0;
			var il = arcLengths.length;
			var targetArcLength; // The targeted u distance value to get

			if (distance) {
				targetArcLength = distance;
			} else {
				targetArcLength = u * arcLengths[il - 1];
			} // binary search for the index with largest value smaller than target u distance


			var low = 0,
					high = il - 1,
					comparison;

			while (low <= high) {
				i = Math.floor(low + (high - low) / 2); // less likely to overflow, though probably not issue here, JS doesn't really have integers, all numbers are floats

				comparison = arcLengths[i] - targetArcLength;

				if (comparison < 0) {
					low = i + 1;
				} else if (comparison > 0) {
					high = i - 1;
				} else {
					high = i;
					break; // DONE
				}
			}

			i = high;

			if (arcLengths[i] === targetArcLength) {
				return i / (il - 1);
			} // we could get finer grain at lengths, or use simple interpolation between two points


			var lengthBefore = arcLengths[i];
			var lengthAfter = arcLengths[i + 1];
			var segmentLength = lengthAfter - lengthBefore; // determine where we are between the 'before' and 'after' points

			var segmentFraction = (targetArcLength - lengthBefore) / segmentLength; // add that fractional amount to t

			var t = (i + segmentFraction) / (il - 1);
			return t;
		},
		// Returns a unit vector tangent at t
		// In case any sub curve does not implement its tangent derivation,
		// 2 points a small delta apart will be used to find its gradient
		// which seems to give a reasonable approximation
		getTangent: function getTangent(t, optionalTarget) {
			var delta = 0.0001;
			var t1 = t - delta;
			var t2 = t + delta; // Capping in case of danger

			if (t1 < 0) t1 = 0;
			if (t2 > 1) t2 = 1;
			var pt1 = this.getPoint(t1);
			var pt2 = this.getPoint(t2);
			var tangent = optionalTarget || (pt1.isVector2 ? new Vector2() : new Vector3());
			tangent.copy(pt2).sub(pt1).normalize();
			return tangent;
		},
		getTangentAt: function getTangentAt(u, optionalTarget) {
			var t = this.getUtoTmapping(u);
			return this.getTangent(t, optionalTarget);
		},
		computeFrenetFrames: function computeFrenetFrames(segments, closed) {
			// see http://www.cs.indiana.edu/pub/techreports/TR425.pdf
			var normal = new Vector3();
			var tangents = [];
			var normals = [];
			var binormals = [];
			var vec = new Vector3();
			var mat = new Matrix4(); // compute the tangent vectors for each segment on the curve

			for (var i = 0; i <= segments; i++) {
				var u = i / segments;
				tangents[i] = this.getTangentAt(u, new Vector3());
				tangents[i].normalize();
			} // select an initial normal vector perpendicular to the first tangent vector,
			// and in the direction of the minimum tangent xyz component


			normals[0] = new Vector3();
			binormals[0] = new Vector3();
			var min = Number.MAX_VALUE;
			var tx = Math.abs(tangents[0].x);
			var ty = Math.abs(tangents[0].y);
			var tz = Math.abs(tangents[0].z);

			if (tx <= min) {
				min = tx;
				normal.set(1, 0, 0);
			}

			if (ty <= min) {
				min = ty;
				normal.set(0, 1, 0);
			}

			if (tz <= min) {
				normal.set(0, 0, 1);
			}

			vec.crossVectors(tangents[0], normal).normalize();
			normals[0].crossVectors(tangents[0], vec);
			binormals[0].crossVectors(tangents[0], normals[0]); // compute the slowly-varying normal and binormal vectors for each segment on the curve

			for (var _i = 1; _i <= segments; _i++) {
				normals[_i] = normals[_i - 1].clone();
				binormals[_i] = binormals[_i - 1].clone();
				vec.crossVectors(tangents[_i - 1], tangents[_i]);

				if (vec.length() > Number.EPSILON) {
					vec.normalize();
					var theta = Math.acos(MathUtils.clamp(tangents[_i - 1].dot(tangents[_i]), -1, 1)); // clamp for floating pt errors

					normals[_i].applyMatrix4(mat.makeRotationAxis(vec, theta));
				}

				binormals[_i].crossVectors(tangents[_i], normals[_i]);
			} // if the curve is closed, postprocess the vectors so the first and last normal vectors are the same


			if (closed === true) {
				var _theta = Math.acos(MathUtils.clamp(normals[0].dot(normals[segments]), -1, 1));

				_theta /= segments;

				if (tangents[0].dot(vec.crossVectors(normals[0], normals[segments])) > 0) {
					_theta = -_theta;
				}

				for (var _i2 = 1; _i2 <= segments; _i2++) {
					// twist a little...
					normals[_i2].applyMatrix4(mat.makeRotationAxis(tangents[_i2], _theta * _i2));

					binormals[_i2].crossVectors(tangents[_i2], normals[_i2]);
				}
			}

			return {
				tangents: tangents,
				normals: normals,
				binormals: binormals
			};
		},
		clone: function clone() {
			return new this.constructor().copy(this);
		},
		copy: function copy(source) {
			this.arcLengthDivisions = source.arcLengthDivisions;
			return this;
		},
		toJSON: function toJSON() {
			var data = {
				metadata: {
					version: 4.5,
					type: 'Curve',
					generator: 'Curve.toJSON'
				}
			};
			data.arcLengthDivisions = this.arcLengthDivisions;
			data.type = this.type;
			return data;
		},
		fromJSON: function fromJSON(json) {
			this.arcLengthDivisions = json.arcLengthDivisions;
			return this;
		}
	});

	function EllipseCurve(aX, aY, xRadius, yRadius, aStartAngle, aEndAngle, aClockwise, aRotation) {
		Curve.call(this);
		this.type = 'EllipseCurve';
		this.aX = aX || 0;
		this.aY = aY || 0;
		this.xRadius = xRadius || 1;
		this.yRadius = yRadius || 1;
		this.aStartAngle = aStartAngle || 0;
		this.aEndAngle = aEndAngle || 2 * Math.PI;
		this.aClockwise = aClockwise || false;
		this.aRotation = aRotation || 0;
	}

	EllipseCurve.prototype = Object.create(Curve.prototype);
	EllipseCurve.prototype.constructor = EllipseCurve;
	EllipseCurve.prototype.isEllipseCurve = true;

	EllipseCurve.prototype.getPoint = function (t, optionalTarget) {
		var point = optionalTarget || new Vector2();
		var twoPi = Math.PI * 2;
		var deltaAngle = this.aEndAngle - this.aStartAngle;
		var samePoints = Math.abs(deltaAngle) < Number.EPSILON; // ensures that deltaAngle is 0 .. 2 PI

		while (deltaAngle < 0) {
			deltaAngle += twoPi;
		}

		while (deltaAngle > twoPi) {
			deltaAngle -= twoPi;
		}

		if (deltaAngle < Number.EPSILON) {
			if (samePoints) {
				deltaAngle = 0;
			} else {
				deltaAngle = twoPi;
			}
		}

		if (this.aClockwise === true && !samePoints) {
			if (deltaAngle === twoPi) {
				deltaAngle = -twoPi;
			} else {
				deltaAngle = deltaAngle - twoPi;
			}
		}

		var angle = this.aStartAngle + t * deltaAngle;
		var x = this.aX + this.xRadius * Math.cos(angle);
		var y = this.aY + this.yRadius * Math.sin(angle);

		if (this.aRotation !== 0) {
			var cos = Math.cos(this.aRotation);
			var sin = Math.sin(this.aRotation);
			var tx = x - this.aX;
			var ty = y - this.aY; // Rotate the point about the center of the ellipse.

			x = tx * cos - ty * sin + this.aX;
			y = tx * sin + ty * cos + this.aY;
		}

		return point.set(x, y);
	};

	EllipseCurve.prototype.copy = function (source) {
		Curve.prototype.copy.call(this, source);
		this.aX = source.aX;
		this.aY = source.aY;
		this.xRadius = source.xRadius;
		this.yRadius = source.yRadius;
		this.aStartAngle = source.aStartAngle;
		this.aEndAngle = source.aEndAngle;
		this.aClockwise = source.aClockwise;
		this.aRotation = source.aRotation;
		return this;
	};

	EllipseCurve.prototype.toJSON = function () {
		var data = Curve.prototype.toJSON.call(this);
		data.aX = this.aX;
		data.aY = this.aY;
		data.xRadius = this.xRadius;
		data.yRadius = this.yRadius;
		data.aStartAngle = this.aStartAngle;
		data.aEndAngle = this.aEndAngle;
		data.aClockwise = this.aClockwise;
		data.aRotation = this.aRotation;
		return data;
	};

	EllipseCurve.prototype.fromJSON = function (json) {
		Curve.prototype.fromJSON.call(this, json);
		this.aX = json.aX;
		this.aY = json.aY;
		this.xRadius = json.xRadius;
		this.yRadius = json.yRadius;
		this.aStartAngle = json.aStartAngle;
		this.aEndAngle = json.aEndAngle;
		this.aClockwise = json.aClockwise;
		this.aRotation = json.aRotation;
		return this;
	};

	function ArcCurve(aX, aY, aRadius, aStartAngle, aEndAngle, aClockwise) {
		EllipseCurve.call(this, aX, aY, aRadius, aRadius, aStartAngle, aEndAngle, aClockwise);
		this.type = 'ArcCurve';
	}

	ArcCurve.prototype = Object.create(EllipseCurve.prototype);
	ArcCurve.prototype.constructor = ArcCurve;
	ArcCurve.prototype.isArcCurve = true;

	/**
	 * Centripetal CatmullRom Curve - which is useful for avoiding
	 * cusps and self-intersections in non-uniform catmull rom curves.
	 * http://www.cemyuksel.com/research/catmullrom_param/catmullrom.pdf
	 *
	 * curve.type accepts centripetal(default), chordal and catmullrom
	 * curve.tension is used for catmullrom which defaults to 0.5
	 */

	/*
	Based on an optimized c++ solution in
	 - http://stackoverflow.com/questions/9489736/catmull-rom-curve-with-no-cusps-and-no-self-intersections/
	 - http://ideone.com/NoEbVM

	This CubicPoly class could be used for reusing some variables and calculations,
	but for three.js curve use, it could be possible inlined and flatten into a single function call
	which can be placed in CurveUtils.
	*/

	function CubicPoly() {
		var c0 = 0,
				c1 = 0,
				c2 = 0,
				c3 = 0;
		/*
		 * Compute coefficients for a cubic polynomial
		 *	 p(s) = c0 + c1*s + c2*s^2 + c3*s^3
		 * such that
		 *	 p(0) = x0, p(1) = x1
		 *	and
		 *	 p'(0) = t0, p'(1) = t1.
		 */

		function init(x0, x1, t0, t1) {
			c0 = x0;
			c1 = t0;
			c2 = -3 * x0 + 3 * x1 - 2 * t0 - t1;
			c3 = 2 * x0 - 2 * x1 + t0 + t1;
		}

		return {
			initCatmullRom: function initCatmullRom(x0, x1, x2, x3, tension) {
				init(x1, x2, tension * (x2 - x0), tension * (x3 - x1));
			},
			initNonuniformCatmullRom: function initNonuniformCatmullRom(x0, x1, x2, x3, dt0, dt1, dt2) {
				// compute tangents when parameterized in [t1,t2]
				var t1 = (x1 - x0) / dt0 - (x2 - x0) / (dt0 + dt1) + (x2 - x1) / dt1;
				var t2 = (x2 - x1) / dt1 - (x3 - x1) / (dt1 + dt2) + (x3 - x2) / dt2; // rescale tangents for parametrization in [0,1]

				t1 *= dt1;
				t2 *= dt1;
				init(x1, x2, t1, t2);
			},
			calc: function calc(t) {
				var t2 = t * t;
				var t3 = t2 * t;
				return c0 + c1 * t + c2 * t2 + c3 * t3;
			}
		};
	} //


	var tmp = new Vector3();
	var px = new CubicPoly(),
			py = new CubicPoly(),
			pz = new CubicPoly();

	function CatmullRomCurve3(points, closed, curveType, tension) {
		Curve.call(this);
		this.type = 'CatmullRomCurve3';
		this.points = points || [];
		this.closed = closed || false;
		this.curveType = curveType || 'centripetal';
		this.tension = tension !== undefined ? tension : 0.5;
	}

	CatmullRomCurve3.prototype = Object.create(Curve.prototype);
	CatmullRomCurve3.prototype.constructor = CatmullRomCurve3;
	CatmullRomCurve3.prototype.isCatmullRomCurve3 = true;

	CatmullRomCurve3.prototype.getPoint = function (t, optionalTarget) {
		var point = optionalTarget || new Vector3();
		var points = this.points;
		var l = points.length;
		var p = (l - (this.closed ? 0 : 1)) * t;
		var intPoint = Math.floor(p);
		var weight = p - intPoint;

		if (this.closed) {
			intPoint += intPoint > 0 ? 0 : (Math.floor(Math.abs(intPoint) / l) + 1) * l;
		} else if (weight === 0 && intPoint === l - 1) {
			intPoint = l - 2;
			weight = 1;
		}

		var p0, p3; // 4 points (p1 & p2 defined below)

		if (this.closed || intPoint > 0) {
			p0 = points[(intPoint - 1) % l];
		} else {
			// extrapolate first point
			tmp.subVectors(points[0], points[1]).add(points[0]);
			p0 = tmp;
		}

		var p1 = points[intPoint % l];
		var p2 = points[(intPoint + 1) % l];

		if (this.closed || intPoint + 2 < l) {
			p3 = points[(intPoint + 2) % l];
		} else {
			// extrapolate last point
			tmp.subVectors(points[l - 1], points[l - 2]).add(points[l - 1]);
			p3 = tmp;
		}

		if (this.curveType === 'centripetal' || this.curveType === 'chordal') {
			// init Centripetal / Chordal Catmull-Rom
			var pow = this.curveType === 'chordal' ? 0.5 : 0.25;
			var dt0 = Math.pow(p0.distanceToSquared(p1), pow);
			var dt1 = Math.pow(p1.distanceToSquared(p2), pow);
			var dt2 = Math.pow(p2.distanceToSquared(p3), pow); // safety check for repeated points

			if (dt1 < 1e-4) dt1 = 1.0;
			if (dt0 < 1e-4) dt0 = dt1;
			if (dt2 < 1e-4) dt2 = dt1;
			px.initNonuniformCatmullRom(p0.x, p1.x, p2.x, p3.x, dt0, dt1, dt2);
			py.initNonuniformCatmullRom(p0.y, p1.y, p2.y, p3.y, dt0, dt1, dt2);
			pz.initNonuniformCatmullRom(p0.z, p1.z, p2.z, p3.z, dt0, dt1, dt2);
		} else if (this.curveType === 'catmullrom') {
			px.initCatmullRom(p0.x, p1.x, p2.x, p3.x, this.tension);
			py.initCatmullRom(p0.y, p1.y, p2.y, p3.y, this.tension);
			pz.initCatmullRom(p0.z, p1.z, p2.z, p3.z, this.tension);
		}

		point.set(px.calc(weight), py.calc(weight), pz.calc(weight));
		return point;
	};

	CatmullRomCurve3.prototype.copy = function (source) {
		Curve.prototype.copy.call(this, source);
		this.points = [];

		for (var i = 0, l = source.points.length; i < l; i++) {
			var point = source.points[i];
			this.points.push(point.clone());
		}

		this.closed = source.closed;
		this.curveType = source.curveType;
		this.tension = source.tension;
		return this;
	};

	CatmullRomCurve3.prototype.toJSON = function () {
		var data = Curve.prototype.toJSON.call(this);
		data.points = [];

		for (var i = 0, l = this.points.length; i < l; i++) {
			var point = this.points[i];
			data.points.push(point.toArray());
		}

		data.closed = this.closed;
		data.curveType = this.curveType;
		data.tension = this.tension;
		return data;
	};

	CatmullRomCurve3.prototype.fromJSON = function (json) {
		Curve.prototype.fromJSON.call(this, json);
		this.points = [];

		for (var i = 0, l = json.points.length; i < l; i++) {
			var point = json.points[i];
			this.points.push(new Vector3().fromArray(point));
		}

		this.closed = json.closed;
		this.curveType = json.curveType;
		this.tension = json.tension;
		return this;
	};

	/**
	 * Bezier Curves formulas obtained from
	 * http://en.wikipedia.org/wiki/Bézier_curve
	 */
	function CatmullRom(t, p0, p1, p2, p3) {
		var v0 = (p2 - p0) * 0.5;
		var v1 = (p3 - p1) * 0.5;
		var t2 = t * t;
		var t3 = t * t2;
		return (2 * p1 - 2 * p2 + v0 + v1) * t3 + (-3 * p1 + 3 * p2 - 2 * v0 - v1) * t2 + v0 * t + p1;
	} //


	function QuadraticBezierP0(t, p) {
		var k = 1 - t;
		return k * k * p;
	}

	function QuadraticBezierP1(t, p) {
		return 2 * (1 - t) * t * p;
	}

	function QuadraticBezierP2(t, p) {
		return t * t * p;
	}

	function QuadraticBezier(t, p0, p1, p2) {
		return QuadraticBezierP0(t, p0) + QuadraticBezierP1(t, p1) + QuadraticBezierP2(t, p2);
	} //


	function CubicBezierP0(t, p) {
		var k = 1 - t;
		return k * k * k * p;
	}

	function CubicBezierP1(t, p) {
		var k = 1 - t;
		return 3 * k * k * t * p;
	}

	function CubicBezierP2(t, p) {
		return 3 * (1 - t) * t * t * p;
	}

	function CubicBezierP3(t, p) {
		return t * t * t * p;
	}

	function CubicBezier(t, p0, p1, p2, p3) {
		return CubicBezierP0(t, p0) + CubicBezierP1(t, p1) + CubicBezierP2(t, p2) + CubicBezierP3(t, p3);
	}

	function CubicBezierCurve(v0, v1, v2, v3) {
		Curve.call(this);
		this.type = 'CubicBezierCurve';
		this.v0 = v0 || new Vector2();
		this.v1 = v1 || new Vector2();
		this.v2 = v2 || new Vector2();
		this.v3 = v3 || new Vector2();
	}

	CubicBezierCurve.prototype = Object.create(Curve.prototype);
	CubicBezierCurve.prototype.constructor = CubicBezierCurve;
	CubicBezierCurve.prototype.isCubicBezierCurve = true;

	CubicBezierCurve.prototype.getPoint = function (t, optionalTarget) {
		var point = optionalTarget || new Vector2();
		var v0 = this.v0,
				v1 = this.v1,
				v2 = this.v2,
				v3 = this.v3;
		point.set(CubicBezier(t, v0.x, v1.x, v2.x, v3.x), CubicBezier(t, v0.y, v1.y, v2.y, v3.y));
		return point;
	};

	CubicBezierCurve.prototype.copy = function (source) {
		Curve.prototype.copy.call(this, source);
		this.v0.copy(source.v0);
		this.v1.copy(source.v1);
		this.v2.copy(source.v2);
		this.v3.copy(source.v3);
		return this;
	};

	CubicBezierCurve.prototype.toJSON = function () {
		var data = Curve.prototype.toJSON.call(this);
		data.v0 = this.v0.toArray();
		data.v1 = this.v1.toArray();
		data.v2 = this.v2.toArray();
		data.v3 = this.v3.toArray();
		return data;
	};

	CubicBezierCurve.prototype.fromJSON = function (json) {
		Curve.prototype.fromJSON.call(this, json);
		this.v0.fromArray(json.v0);
		this.v1.fromArray(json.v1);
		this.v2.fromArray(json.v2);
		this.v3.fromArray(json.v3);
		return this;
	};

	function CubicBezierCurve3(v0, v1, v2, v3) {
		Curve.call(this);
		this.type = 'CubicBezierCurve3';
		this.v0 = v0 || new Vector3();
		this.v1 = v1 || new Vector3();
		this.v2 = v2 || new Vector3();
		this.v3 = v3 || new Vector3();
	}

	CubicBezierCurve3.prototype = Object.create(Curve.prototype);
	CubicBezierCurve3.prototype.constructor = CubicBezierCurve3;
	CubicBezierCurve3.prototype.isCubicBezierCurve3 = true;

	CubicBezierCurve3.prototype.getPoint = function (t, optionalTarget) {
		var point = optionalTarget || new Vector3();
		var v0 = this.v0,
				v1 = this.v1,
				v2 = this.v2,
				v3 = this.v3;
		point.set(CubicBezier(t, v0.x, v1.x, v2.x, v3.x), CubicBezier(t, v0.y, v1.y, v2.y, v3.y), CubicBezier(t, v0.z, v1.z, v2.z, v3.z));
		return point;
	};

	CubicBezierCurve3.prototype.copy = function (source) {
		Curve.prototype.copy.call(this, source);
		this.v0.copy(source.v0);
		this.v1.copy(source.v1);
		this.v2.copy(source.v2);
		this.v3.copy(source.v3);
		return this;
	};

	CubicBezierCurve3.prototype.toJSON = function () {
		var data = Curve.prototype.toJSON.call(this);
		data.v0 = this.v0.toArray();
		data.v1 = this.v1.toArray();
		data.v2 = this.v2.toArray();
		data.v3 = this.v3.toArray();
		return data;
	};

	CubicBezierCurve3.prototype.fromJSON = function (json) {
		Curve.prototype.fromJSON.call(this, json);
		this.v0.fromArray(json.v0);
		this.v1.fromArray(json.v1);
		this.v2.fromArray(json.v2);
		this.v3.fromArray(json.v3);
		return this;
	};

	function LineCurve(v1, v2) {
		Curve.call(this);
		this.type = 'LineCurve';
		this.v1 = v1 || new Vector2();
		this.v2 = v2 || new Vector2();
	}

	LineCurve.prototype = Object.create(Curve.prototype);
	LineCurve.prototype.constructor = LineCurve;
	LineCurve.prototype.isLineCurve = true;

	LineCurve.prototype.getPoint = function (t, optionalTarget) {
		var point = optionalTarget || new Vector2();

		if (t === 1) {
			point.copy(this.v2);
		} else {
			point.copy(this.v2).sub(this.v1);
			point.multiplyScalar(t).add(this.v1);
		}

		return point;
	}; // Line curve is linear, so we can overwrite default getPointAt


	LineCurve.prototype.getPointAt = function (u, optionalTarget) {
		return this.getPoint(u, optionalTarget);
	};

	LineCurve.prototype.getTangent = function (t, optionalTarget) {
		var tangent = optionalTarget || new Vector2();
		tangent.copy(this.v2).sub(this.v1).normalize();
		return tangent;
	};

	LineCurve.prototype.copy = function (source) {
		Curve.prototype.copy.call(this, source);
		this.v1.copy(source.v1);
		this.v2.copy(source.v2);
		return this;
	};

	LineCurve.prototype.toJSON = function () {
		var data = Curve.prototype.toJSON.call(this);
		data.v1 = this.v1.toArray();
		data.v2 = this.v2.toArray();
		return data;
	};

	LineCurve.prototype.fromJSON = function (json) {
		Curve.prototype.fromJSON.call(this, json);
		this.v1.fromArray(json.v1);
		this.v2.fromArray(json.v2);
		return this;
	};

	function LineCurve3(v1, v2) {
		Curve.call(this);
		this.type = 'LineCurve3';
		this.v1 = v1 || new Vector3();
		this.v2 = v2 || new Vector3();
	}

	LineCurve3.prototype = Object.create(Curve.prototype);
	LineCurve3.prototype.constructor = LineCurve3;
	LineCurve3.prototype.isLineCurve3 = true;

	LineCurve3.prototype.getPoint = function (t, optionalTarget) {
		var point = optionalTarget || new Vector3();

		if (t === 1) {
			point.copy(this.v2);
		} else {
			point.copy(this.v2).sub(this.v1);
			point.multiplyScalar(t).add(this.v1);
		}

		return point;
	}; // Line curve is linear, so we can overwrite default getPointAt


	LineCurve3.prototype.getPointAt = function (u, optionalTarget) {
		return this.getPoint(u, optionalTarget);
	};

	LineCurve3.prototype.copy = function (source) {
		Curve.prototype.copy.call(this, source);
		this.v1.copy(source.v1);
		this.v2.copy(source.v2);
		return this;
	};

	LineCurve3.prototype.toJSON = function () {
		var data = Curve.prototype.toJSON.call(this);
		data.v1 = this.v1.toArray();
		data.v2 = this.v2.toArray();
		return data;
	};

	LineCurve3.prototype.fromJSON = function (json) {
		Curve.prototype.fromJSON.call(this, json);
		this.v1.fromArray(json.v1);
		this.v2.fromArray(json.v2);
		return this;
	};

	function QuadraticBezierCurve(v0, v1, v2) {
		Curve.call(this);
		this.type = 'QuadraticBezierCurve';
		this.v0 = v0 || new Vector2();
		this.v1 = v1 || new Vector2();
		this.v2 = v2 || new Vector2();
	}

	QuadraticBezierCurve.prototype = Object.create(Curve.prototype);
	QuadraticBezierCurve.prototype.constructor = QuadraticBezierCurve;
	QuadraticBezierCurve.prototype.isQuadraticBezierCurve = true;

	QuadraticBezierCurve.prototype.getPoint = function (t, optionalTarget) {
		var point = optionalTarget || new Vector2();
		var v0 = this.v0,
				v1 = this.v1,
				v2 = this.v2;
		point.set(QuadraticBezier(t, v0.x, v1.x, v2.x), QuadraticBezier(t, v0.y, v1.y, v2.y));
		return point;
	};

	QuadraticBezierCurve.prototype.copy = function (source) {
		Curve.prototype.copy.call(this, source);
		this.v0.copy(source.v0);
		this.v1.copy(source.v1);
		this.v2.copy(source.v2);
		return this;
	};

	QuadraticBezierCurve.prototype.toJSON = function () {
		var data = Curve.prototype.toJSON.call(this);
		data.v0 = this.v0.toArray();
		data.v1 = this.v1.toArray();
		data.v2 = this.v2.toArray();
		return data;
	};

	QuadraticBezierCurve.prototype.fromJSON = function (json) {
		Curve.prototype.fromJSON.call(this, json);
		this.v0.fromArray(json.v0);
		this.v1.fromArray(json.v1);
		this.v2.fromArray(json.v2);
		return this;
	};

	function QuadraticBezierCurve3(v0, v1, v2) {
		Curve.call(this);
		this.type = 'QuadraticBezierCurve3';
		this.v0 = v0 || new Vector3();
		this.v1 = v1 || new Vector3();
		this.v2 = v2 || new Vector3();
	}

	QuadraticBezierCurve3.prototype = Object.create(Curve.prototype);
	QuadraticBezierCurve3.prototype.constructor = QuadraticBezierCurve3;
	QuadraticBezierCurve3.prototype.isQuadraticBezierCurve3 = true;

	QuadraticBezierCurve3.prototype.getPoint = function (t, optionalTarget) {
		var point = optionalTarget || new Vector3();
		var v0 = this.v0,
				v1 = this.v1,
				v2 = this.v2;
		point.set(QuadraticBezier(t, v0.x, v1.x, v2.x), QuadraticBezier(t, v0.y, v1.y, v2.y), QuadraticBezier(t, v0.z, v1.z, v2.z));
		return point;
	};

	QuadraticBezierCurve3.prototype.copy = function (source) {
		Curve.prototype.copy.call(this, source);
		this.v0.copy(source.v0);
		this.v1.copy(source.v1);
		this.v2.copy(source.v2);
		return this;
	};

	QuadraticBezierCurve3.prototype.toJSON = function () {
		var data = Curve.prototype.toJSON.call(this);
		data.v0 = this.v0.toArray();
		data.v1 = this.v1.toArray();
		data.v2 = this.v2.toArray();
		return data;
	};

	QuadraticBezierCurve3.prototype.fromJSON = function (json) {
		Curve.prototype.fromJSON.call(this, json);
		this.v0.fromArray(json.v0);
		this.v1.fromArray(json.v1);
		this.v2.fromArray(json.v2);
		return this;
	};

	function SplineCurve(points) {
		Curve.call(this);
		this.type = 'SplineCurve';
		this.points = points || [];
	}

	SplineCurve.prototype = Object.create(Curve.prototype);
	SplineCurve.prototype.constructor = SplineCurve;
	SplineCurve.prototype.isSplineCurve = true;

	SplineCurve.prototype.getPoint = function (t, optionalTarget) {
		var point = optionalTarget || new Vector2();
		var points = this.points;
		var p = (points.length - 1) * t;
		var intPoint = Math.floor(p);
		var weight = p - intPoint;
		var p0 = points[intPoint === 0 ? intPoint : intPoint - 1];
		var p1 = points[intPoint];
		var p2 = points[intPoint > points.length - 2 ? points.length - 1 : intPoint + 1];
		var p3 = points[intPoint > points.length - 3 ? points.length - 1 : intPoint + 2];
		point.set(CatmullRom(weight, p0.x, p1.x, p2.x, p3.x), CatmullRom(weight, p0.y, p1.y, p2.y, p3.y));
		return point;
	};

	SplineCurve.prototype.copy = function (source) {
		Curve.prototype.copy.call(this, source);
		this.points = [];

		for (var i = 0, l = source.points.length; i < l; i++) {
			var point = source.points[i];
			this.points.push(point.clone());
		}

		return this;
	};

	SplineCurve.prototype.toJSON = function () {
		var data = Curve.prototype.toJSON.call(this);
		data.points = [];

		for (var i = 0, l = this.points.length; i < l; i++) {
			var point = this.points[i];
			data.points.push(point.toArray());
		}

		return data;
	};

	SplineCurve.prototype.fromJSON = function (json) {
		Curve.prototype.fromJSON.call(this, json);
		this.points = [];

		for (var i = 0, l = json.points.length; i < l; i++) {
			var point = json.points[i];
			this.points.push(new Vector2().fromArray(point));
		}

		return this;
	};

	var Curves = /*#__PURE__*/Object.freeze({
		__proto__: null,
		ArcCurve: ArcCurve,
		CatmullRomCurve3: CatmullRomCurve3,
		CubicBezierCurve: CubicBezierCurve,
		CubicBezierCurve3: CubicBezierCurve3,
		EllipseCurve: EllipseCurve,
		LineCurve: LineCurve,
		LineCurve3: LineCurve3,
		QuadraticBezierCurve: QuadraticBezierCurve,
		QuadraticBezierCurve3: QuadraticBezierCurve3,
		SplineCurve: SplineCurve
	});

	/**************************************************************
	 *	Curved Path - a curve path is simply a array of connected
	 *	curves, but retains the api of a curve
	 **************************************************************/

	function CurvePath() {
		Curve.call(this);
		this.type = 'CurvePath';
		this.curves = [];
		this.autoClose = false; // Automatically closes the path
	}

	CurvePath.prototype = Object.assign(Object.create(Curve.prototype), {
		constructor: CurvePath,
		add: function add(curve) {
			this.curves.push(curve);
		},
		closePath: function closePath() {
			// Add a line curve if start and end of lines are not connected
			var startPoint = this.curves[0].getPoint(0);
			var endPoint = this.curves[this.curves.length - 1].getPoint(1);

			if (!startPoint.equals(endPoint)) {
				this.curves.push(new LineCurve(endPoint, startPoint));
			}
		},
		// To get accurate point with reference to
		// entire path distance at time t,
		// following has to be done:
		// 1. Length of each sub path have to be known
		// 2. Locate and identify type of curve
		// 3. Get t for the curve
		// 4. Return curve.getPointAt(t')
		getPoint: function getPoint(t) {
			var d = t * this.getLength();
			var curveLengths = this.getCurveLengths();
			var i = 0; // To think about boundaries points.

			while (i < curveLengths.length) {
				if (curveLengths[i] >= d) {
					var diff = curveLengths[i] - d;
					var curve = this.curves[i];
					var segmentLength = curve.getLength();
					var u = segmentLength === 0 ? 0 : 1 - diff / segmentLength;
					return curve.getPointAt(u);
				}

				i++;
			}

			return null; // loop where sum != 0, sum > d , sum+1 <d
		},
		// We cannot use the default THREE.Curve getPoint() with getLength() because in
		// THREE.Curve, getLength() depends on getPoint() but in THREE.CurvePath
		// getPoint() depends on getLength
		getLength: function getLength() {
			var lens = this.getCurveLengths();
			return lens[lens.length - 1];
		},
		// cacheLengths must be recalculated.
		updateArcLengths: function updateArcLengths() {
			this.needsUpdate = true;
			this.cacheLengths = null;
			this.getCurveLengths();
		},
		// Compute lengths and cache them
		// We cannot overwrite getLengths() because UtoT mapping uses it.
		getCurveLengths: function getCurveLengths() {
			// We use cache values if curves and cache array are same length
			if (this.cacheLengths && this.cacheLengths.length === this.curves.length) {
				return this.cacheLengths;
			} // Get length of sub-curve
			// Push sums into cached array


			var lengths = [];
			var sums = 0;

			for (var i = 0, l = this.curves.length; i < l; i++) {
				sums += this.curves[i].getLength();
				lengths.push(sums);
			}

			this.cacheLengths = lengths;
			return lengths;
		},
		getSpacedPoints: function getSpacedPoints(divisions) {
			if (divisions === undefined) divisions = 40;
			var points = [];

			for (var i = 0; i <= divisions; i++) {
				points.push(this.getPoint(i / divisions));
			}

			if (this.autoClose) {
				points.push(points[0]);
			}

			return points;
		},
		getPoints: function getPoints(divisions) {
			divisions = divisions || 12;
			var points = [];
			var last;

			for (var i = 0, curves = this.curves; i < curves.length; i++) {
				var curve = curves[i];
				var resolution = curve && curve.isEllipseCurve ? divisions * 2 : curve && (curve.isLineCurve || curve.isLineCurve3) ? 1 : curve && curve.isSplineCurve ? divisions * curve.points.length : divisions;
				var pts = curve.getPoints(resolution);

				for (var j = 0; j < pts.length; j++) {
					var point = pts[j];
					if (last && last.equals(point)) continue; // ensures no consecutive points are duplicates

					points.push(point);
					last = point;
				}
			}

			if (this.autoClose && points.length > 1 && !points[points.length - 1].equals(points[0])) {
				points.push(points[0]);
			}

			return points;
		},
		copy: function copy(source) {
			Curve.prototype.copy.call(this, source);
			this.curves = [];

			for (var i = 0, l = source.curves.length; i < l; i++) {
				var curve = source.curves[i];
				this.curves.push(curve.clone());
			}

			this.autoClose = source.autoClose;
			return this;
		},
		toJSON: function toJSON() {
			var data = Curve.prototype.toJSON.call(this);
			data.autoClose = this.autoClose;
			data.curves = [];

			for (var i = 0, l = this.curves.length; i < l; i++) {
				var curve = this.curves[i];
				data.curves.push(curve.toJSON());
			}

			return data;
		},
		fromJSON: function fromJSON(json) {
			Curve.prototype.fromJSON.call(this, json);
			this.autoClose = json.autoClose;
			this.curves = [];

			for (var i = 0, l = json.curves.length; i < l; i++) {
				var curve = json.curves[i];
				this.curves.push(new Curves[curve.type]().fromJSON(curve));
			}

			return this;
		}
	});

	function Path(points) {
		CurvePath.call(this);
		this.type = 'Path';
		this.currentPoint = new Vector2();

		if (points) {
			this.setFromPoints(points);
		}
	}

	Path.prototype = Object.assign(Object.create(CurvePath.prototype), {
		constructor: Path,
		setFromPoints: function setFromPoints(points) {
			this.moveTo(points[0].x, points[0].y);

			for (var i = 1, l = points.length; i < l; i++) {
				this.lineTo(points[i].x, points[i].y);
			}

			return this;
		},
		moveTo: function moveTo(x, y) {
			this.currentPoint.set(x, y); // TODO consider referencing vectors instead of copying?

			return this;
		},
		lineTo: function lineTo(x, y) {
			var curve = new LineCurve(this.currentPoint.clone(), new Vector2(x, y));
			this.curves.push(curve);
			this.currentPoint.set(x, y);
			return this;
		},
		quadraticCurveTo: function quadraticCurveTo(aCPx, aCPy, aX, aY) {
			var curve = new QuadraticBezierCurve(this.currentPoint.clone(), new Vector2(aCPx, aCPy), new Vector2(aX, aY));
			this.curves.push(curve);
			this.currentPoint.set(aX, aY);
			return this;
		},
		bezierCurveTo: function bezierCurveTo(aCP1x, aCP1y, aCP2x, aCP2y, aX, aY) {
			var curve = new CubicBezierCurve(this.currentPoint.clone(), new Vector2(aCP1x, aCP1y), new Vector2(aCP2x, aCP2y), new Vector2(aX, aY));
			this.curves.push(curve);
			this.currentPoint.set(aX, aY);
			return this;
		},
		splineThru: function splineThru(pts
		/*Array of Vector*/
		) {
			var npts = [this.currentPoint.clone()].concat(pts);
			var curve = new SplineCurve(npts);
			this.curves.push(curve);
			this.currentPoint.copy(pts[pts.length - 1]);
			return this;
		},
		arc: function arc(aX, aY, aRadius, aStartAngle, aEndAngle, aClockwise) {
			var x0 = this.currentPoint.x;
			var y0 = this.currentPoint.y;
			this.absarc(aX + x0, aY + y0, aRadius, aStartAngle, aEndAngle, aClockwise);
			return this;
		},
		absarc: function absarc(aX, aY, aRadius, aStartAngle, aEndAngle, aClockwise) {
			this.absellipse(aX, aY, aRadius, aRadius, aStartAngle, aEndAngle, aClockwise);
			return this;
		},
		ellipse: function ellipse(aX, aY, xRadius, yRadius, aStartAngle, aEndAngle, aClockwise, aRotation) {
			var x0 = this.currentPoint.x;
			var y0 = this.currentPoint.y;
			this.absellipse(aX + x0, aY + y0, xRadius, yRadius, aStartAngle, aEndAngle, aClockwise, aRotation);
			return this;
		},
		absellipse: function absellipse(aX, aY, xRadius, yRadius, aStartAngle, aEndAngle, aClockwise, aRotation) {
			var curve = new EllipseCurve(aX, aY, xRadius, yRadius, aStartAngle, aEndAngle, aClockwise, aRotation);

			if (this.curves.length > 0) {
				// if a previous curve is present, attempt to join
				var firstPoint = curve.getPoint(0);

				if (!firstPoint.equals(this.currentPoint)) {
					this.lineTo(firstPoint.x, firstPoint.y);
				}
			}

			this.curves.push(curve);
			var lastPoint = curve.getPoint(1);
			this.currentPoint.copy(lastPoint);
			return this;
		},
		copy: function copy(source) {
			CurvePath.prototype.copy.call(this, source);
			this.currentPoint.copy(source.currentPoint);
			return this;
		},
		toJSON: function toJSON() {
			var data = CurvePath.prototype.toJSON.call(this);
			data.currentPoint = this.currentPoint.toArray();
			return data;
		},
		fromJSON: function fromJSON(json) {
			CurvePath.prototype.fromJSON.call(this, json);
			this.currentPoint.fromArray(json.currentPoint);
			return this;
		}
	});

	function Shape(points) {
		Path.call(this, points);
		this.uuid = MathUtils.generateUUID();
		this.type = 'Shape';
		this.holes = [];
	}

	Shape.prototype = Object.assign(Object.create(Path.prototype), {
		constructor: Shape,
		getPointsHoles: function getPointsHoles(divisions) {
			var holesPts = [];

			for (var i = 0, l = this.holes.length; i < l; i++) {
				holesPts[i] = this.holes[i].getPoints(divisions);
			}

			return holesPts;
		},
		// get points of shape and holes (keypoints based on segments parameter)
		extractPoints: function extractPoints(divisions) {
			return {
				shape: this.getPoints(divisions),
				holes: this.getPointsHoles(divisions)
			};
		},
		copy: function copy(source) {
			Path.prototype.copy.call(this, source);
			this.holes = [];

			for (var i = 0, l = source.holes.length; i < l; i++) {
				var hole = source.holes[i];
				this.holes.push(hole.clone());
			}

			return this;
		},
		toJSON: function toJSON() {
			var data = Path.prototype.toJSON.call(this);
			data.uuid = this.uuid;
			data.holes = [];

			for (var i = 0, l = this.holes.length; i < l; i++) {
				var hole = this.holes[i];
				data.holes.push(hole.toJSON());
			}

			return data;
		},
		fromJSON: function fromJSON(json) {
			Path.prototype.fromJSON.call(this, json);
			this.uuid = json.uuid;
			this.holes = [];

			for (var i = 0, l = json.holes.length; i < l; i++) {
				var hole = json.holes[i];
				this.holes.push(new Path().fromJSON(hole));
			}

			return this;
		}
	});

	var _colorKeywords = {
		'aliceblue': 0xF0F8FF,
		'antiquewhite': 0xFAEBD7,
		'aqua': 0x00FFFF,
		'aquamarine': 0x7FFFD4,
		'azure': 0xF0FFFF,
		'beige': 0xF5F5DC,
		'bisque': 0xFFE4C4,
		'black': 0x000000,
		'blanchedalmond': 0xFFEBCD,
		'blue': 0x0000FF,
		'blueviolet': 0x8A2BE2,
		'brown': 0xA52A2A,
		'burlywood': 0xDEB887,
		'cadetblue': 0x5F9EA0,
		'chartreuse': 0x7FFF00,
		'chocolate': 0xD2691E,
		'coral': 0xFF7F50,
		'cornflowerblue': 0x6495ED,
		'cornsilk': 0xFFF8DC,
		'crimson': 0xDC143C,
		'cyan': 0x00FFFF,
		'darkblue': 0x00008B,
		'darkcyan': 0x008B8B,
		'darkgoldenrod': 0xB8860B,
		'darkgray': 0xA9A9A9,
		'darkgreen': 0x006400,
		'darkgrey': 0xA9A9A9,
		'darkkhaki': 0xBDB76B,
		'darkmagenta': 0x8B008B,
		'darkolivegreen': 0x556B2F,
		'darkorange': 0xFF8C00,
		'darkorchid': 0x9932CC,
		'darkred': 0x8B0000,
		'darksalmon': 0xE9967A,
		'darkseagreen': 0x8FBC8F,
		'darkslateblue': 0x483D8B,
		'darkslategray': 0x2F4F4F,
		'darkslategrey': 0x2F4F4F,
		'darkturquoise': 0x00CED1,
		'darkviolet': 0x9400D3,
		'deeppink': 0xFF1493,
		'deepskyblue': 0x00BFFF,
		'dimgray': 0x696969,
		'dimgrey': 0x696969,
		'dodgerblue': 0x1E90FF,
		'firebrick': 0xB22222,
		'floralwhite': 0xFFFAF0,
		'forestgreen': 0x228B22,
		'fuchsia': 0xFF00FF,
		'gainsboro': 0xDCDCDC,
		'ghostwhite': 0xF8F8FF,
		'gold': 0xFFD700,
		'goldenrod': 0xDAA520,
		'gray': 0x808080,
		'green': 0x008000,
		'greenyellow': 0xADFF2F,
		'grey': 0x808080,
		'honeydew': 0xF0FFF0,
		'hotpink': 0xFF69B4,
		'indianred': 0xCD5C5C,
		'indigo': 0x4B0082,
		'ivory': 0xFFFFF0,
		'khaki': 0xF0E68C,
		'lavender': 0xE6E6FA,
		'lavenderblush': 0xFFF0F5,
		'lawngreen': 0x7CFC00,
		'lemonchiffon': 0xFFFACD,
		'lightblue': 0xADD8E6,
		'lightcoral': 0xF08080,
		'lightcyan': 0xE0FFFF,
		'lightgoldenrodyellow': 0xFAFAD2,
		'lightgray': 0xD3D3D3,
		'lightgreen': 0x90EE90,
		'lightgrey': 0xD3D3D3,
		'lightpink': 0xFFB6C1,
		'lightsalmon': 0xFFA07A,
		'lightseagreen': 0x20B2AA,
		'lightskyblue': 0x87CEFA,
		'lightslategray': 0x778899,
		'lightslategrey': 0x778899,
		'lightsteelblue': 0xB0C4DE,
		'lightyellow': 0xFFFFE0,
		'lime': 0x00FF00,
		'limegreen': 0x32CD32,
		'linen': 0xFAF0E6,
		'magenta': 0xFF00FF,
		'maroon': 0x800000,
		'mediumaquamarine': 0x66CDAA,
		'mediumblue': 0x0000CD,
		'mediumorchid': 0xBA55D3,
		'mediumpurple': 0x9370DB,
		'mediumseagreen': 0x3CB371,
		'mediumslateblue': 0x7B68EE,
		'mediumspringgreen': 0x00FA9A,
		'mediumturquoise': 0x48D1CC,
		'mediumvioletred': 0xC71585,
		'midnightblue': 0x191970,
		'mintcream': 0xF5FFFA,
		'mistyrose': 0xFFE4E1,
		'moccasin': 0xFFE4B5,
		'navajowhite': 0xFFDEAD,
		'navy': 0x000080,
		'oldlace': 0xFDF5E6,
		'olive': 0x808000,
		'olivedrab': 0x6B8E23,
		'orange': 0xFFA500,
		'orangered': 0xFF4500,
		'orchid': 0xDA70D6,
		'palegoldenrod': 0xEEE8AA,
		'palegreen': 0x98FB98,
		'paleturquoise': 0xAFEEEE,
		'palevioletred': 0xDB7093,
		'papayawhip': 0xFFEFD5,
		'peachpuff': 0xFFDAB9,
		'peru': 0xCD853F,
		'pink': 0xFFC0CB,
		'plum': 0xDDA0DD,
		'powderblue': 0xB0E0E6,
		'purple': 0x800080,
		'rebeccapurple': 0x663399,
		'red': 0xFF0000,
		'rosybrown': 0xBC8F8F,
		'royalblue': 0x4169E1,
		'saddlebrown': 0x8B4513,
		'salmon': 0xFA8072,
		'sandybrown': 0xF4A460,
		'seagreen': 0x2E8B57,
		'seashell': 0xFFF5EE,
		'sienna': 0xA0522D,
		'silver': 0xC0C0C0,
		'skyblue': 0x87CEEB,
		'slateblue': 0x6A5ACD,
		'slategray': 0x708090,
		'slategrey': 0x708090,
		'snow': 0xFFFAFA,
		'springgreen': 0x00FF7F,
		'steelblue': 0x4682B4,
		'tan': 0xD2B48C,
		'teal': 0x008080,
		'thistle': 0xD8BFD8,
		'tomato': 0xFF6347,
		'turquoise': 0x40E0D0,
		'violet': 0xEE82EE,
		'wheat': 0xF5DEB3,
		'white': 0xFFFFFF,
		'whitesmoke': 0xF5F5F5,
		'yellow': 0xFFFF00,
		'yellowgreen': 0x9ACD32
	};
	var _hslA = {
		h: 0,
		s: 0,
		l: 0
	};
	var _hslB = {
		h: 0,
		s: 0,
		l: 0
	};

	function hue2rgb(p, q, t) {
		if (t < 0) t += 1;
		if (t > 1) t -= 1;
		if (t < 1 / 6) return p + (q - p) * 6 * t;
		if (t < 1 / 2) return q;
		if (t < 2 / 3) return p + (q - p) * 6 * (2 / 3 - t);
		return p;
	}

	function SRGBToLinear(c) {
		return c < 0.04045 ? c * 0.0773993808 : Math.pow(c * 0.9478672986 + 0.0521327014, 2.4);
	}

	function LinearToSRGB(c) {
		return c < 0.0031308 ? c * 12.92 : 1.055 * Math.pow(c, 0.41666) - 0.055;
	}

	var Color = /*#__PURE__*/function () {
		function Color(r, g, b) {
			Object.defineProperty(this, 'isColor', {
				value: true
			});

			if (g === undefined && b === undefined) {
				// r is THREE.Color, hex or string
				return this.set(r);
			}

			return this.setRGB(r, g, b);
		}

		var _proto = Color.prototype;

		_proto.set = function set(value) {
			if (value && value.isColor) {
				this.copy(value);
			} else if (typeof value === 'number') {
				this.setHex(value);
			} else if (typeof value === 'string') {
				this.setStyle(value);
			}

			return this;
		};

		_proto.setScalar = function setScalar(scalar) {
			this.r = scalar;
			this.g = scalar;
			this.b = scalar;
			return this;
		};

		_proto.setHex = function setHex(hex) {
			hex = Math.floor(hex);
			this.r = (hex >> 16 & 255) / 255;
			this.g = (hex >> 8 & 255) / 255;
			this.b = (hex & 255) / 255;
			return this;
		};

		_proto.setRGB = function setRGB(r, g, b) {
			this.r = r;
			this.g = g;
			this.b = b;
			return this;
		};

		_proto.setHSL = function setHSL(h, s, l) {
			// h,s,l ranges are in 0.0 - 1.0
			h = MathUtils.euclideanModulo(h, 1);
			s = MathUtils.clamp(s, 0, 1);
			l = MathUtils.clamp(l, 0, 1);

			if (s === 0) {
				this.r = this.g = this.b = l;
			} else {
				var p = l <= 0.5 ? l * (1 + s) : l + s - l * s;
				var q = 2 * l - p;
				this.r = hue2rgb(q, p, h + 1 / 3);
				this.g = hue2rgb(q, p, h);
				this.b = hue2rgb(q, p, h - 1 / 3);
			}

			return this;
		};

		_proto.setStyle = function setStyle(style) {
			function handleAlpha(string) {
				if (string === undefined) return;

				if (parseFloat(string) < 1) {
					console.warn('THREE.Color: Alpha component of ' + style + ' will be ignored.');
				}
			}

			var m;

			if (m = /^((?:rgb|hsl)a?)\(\s*([^\)]*)\)/.exec(style)) {
				// rgb / hsl
				var color;
				var name = m[1];
				var components = m[2];

				switch (name) {
					case 'rgb':
					case 'rgba':
						if (color = /^(\d+)\s*,\s*(\d+)\s*,\s*(\d+)\s*(,\s*([0-9]*\.?[0-9]+)\s*)?$/.exec(components)) {
							// rgb(255,0,0) rgba(255,0,0,0.5)
							this.r = Math.min(255, parseInt(color[1], 10)) / 255;
							this.g = Math.min(255, parseInt(color[2], 10)) / 255;
							this.b = Math.min(255, parseInt(color[3], 10)) / 255;
							handleAlpha(color[5]);
							return this;
						}

						if (color = /^(\d+)\%\s*,\s*(\d+)\%\s*,\s*(\d+)\%\s*(,\s*([0-9]*\.?[0-9]+)\s*)?$/.exec(components)) {
							// rgb(100%,0%,0%) rgba(100%,0%,0%,0.5)
							this.r = Math.min(100, parseInt(color[1], 10)) / 100;
							this.g = Math.min(100, parseInt(color[2], 10)) / 100;
							this.b = Math.min(100, parseInt(color[3], 10)) / 100;
							handleAlpha(color[5]);
							return this;
						}

						break;

					case 'hsl':
					case 'hsla':
						if (color = /^([0-9]*\.?[0-9]+)\s*,\s*(\d+)\%\s*,\s*(\d+)\%\s*(,\s*([0-9]*\.?[0-9]+)\s*)?$/.exec(components)) {
							// hsl(120,50%,50%) hsla(120,50%,50%,0.5)
							var h = parseFloat(color[1]) / 360;
							var s = parseInt(color[2], 10) / 100;
							var l = parseInt(color[3], 10) / 100;
							handleAlpha(color[5]);
							return this.setHSL(h, s, l);
						}

						break;
				}
			} else if (m = /^\#([A-Fa-f0-9]+)$/.exec(style)) {
				// hex color
				var hex = m[1];
				var size = hex.length;

				if (size === 3) {
					// #ff0
					this.r = parseInt(hex.charAt(0) + hex.charAt(0), 16) / 255;
					this.g = parseInt(hex.charAt(1) + hex.charAt(1), 16) / 255;
					this.b = parseInt(hex.charAt(2) + hex.charAt(2), 16) / 255;
					return this;
				} else if (size === 6) {
					// #ff0000
					this.r = parseInt(hex.charAt(0) + hex.charAt(1), 16) / 255;
					this.g = parseInt(hex.charAt(2) + hex.charAt(3), 16) / 255;
					this.b = parseInt(hex.charAt(4) + hex.charAt(5), 16) / 255;
					return this;
				}
			}

			if (style && style.length > 0) {
				return this.setColorName(style);
			}

			return this;
		};

		_proto.setColorName = function setColorName(style) {
			// color keywords
			var hex = _colorKeywords[style];

			if (hex !== undefined) {
				// red
				this.setHex(hex);
			} else {
				// unknown color
				console.warn('THREE.Color: Unknown color ' + style);
			}

			return this;
		};

		_proto.clone = function clone() {
			return new this.constructor(this.r, this.g, this.b);
		};

		_proto.copy = function copy(color) {
			this.r = color.r;
			this.g = color.g;
			this.b = color.b;
			return this;
		};

		_proto.copyGammaToLinear = function copyGammaToLinear(color, gammaFactor) {
			if (gammaFactor === undefined) gammaFactor = 2.0;
			this.r = Math.pow(color.r, gammaFactor);
			this.g = Math.pow(color.g, gammaFactor);
			this.b = Math.pow(color.b, gammaFactor);
			return this;
		};

		_proto.copyLinearToGamma = function copyLinearToGamma(color, gammaFactor) {
			if (gammaFactor === undefined) gammaFactor = 2.0;
			var safeInverse = gammaFactor > 0 ? 1.0 / gammaFactor : 1.0;
			this.r = Math.pow(color.r, safeInverse);
			this.g = Math.pow(color.g, safeInverse);
			this.b = Math.pow(color.b, safeInverse);
			return this;
		};

		_proto.convertGammaToLinear = function convertGammaToLinear(gammaFactor) {
			this.copyGammaToLinear(this, gammaFactor);
			return this;
		};

		_proto.convertLinearToGamma = function convertLinearToGamma(gammaFactor) {
			this.copyLinearToGamma(this, gammaFactor);
			return this;
		};

		_proto.copySRGBToLinear = function copySRGBToLinear(color) {
			this.r = SRGBToLinear(color.r);
			this.g = SRGBToLinear(color.g);
			this.b = SRGBToLinear(color.b);
			return this;
		};

		_proto.copyLinearToSRGB = function copyLinearToSRGB(color) {
			this.r = LinearToSRGB(color.r);
			this.g = LinearToSRGB(color.g);
			this.b = LinearToSRGB(color.b);
			return this;
		};

		_proto.convertSRGBToLinear = function convertSRGBToLinear() {
			this.copySRGBToLinear(this);
			return this;
		};

		_proto.convertLinearToSRGB = function convertLinearToSRGB() {
			this.copyLinearToSRGB(this);
			return this;
		};

		_proto.getHex = function getHex() {
			return this.r * 255 << 16 ^ this.g * 255 << 8 ^ this.b * 255 << 0;
		};

		_proto.getHexString = function getHexString() {
			return ('000000' + this.getHex().toString(16)).slice(-6);
		};

		_proto.getHSL = function getHSL(target) {
			// h,s,l ranges are in 0.0 - 1.0
			if (target === undefined) {
				console.warn('THREE.Color: .getHSL() target is now required');
				target = {
					h: 0,
					s: 0,
					l: 0
				};
			}

			var r = this.r,
					g = this.g,
					b = this.b;
			var max = Math.max(r, g, b);
			var min = Math.min(r, g, b);
			var hue, saturation;
			var lightness = (min + max) / 2.0;

			if (min === max) {
				hue = 0;
				saturation = 0;
			} else {
				var delta = max - min;
				saturation = lightness <= 0.5 ? delta / (max + min) : delta / (2 - max - min);

				switch (max) {
					case r:
						hue = (g - b) / delta + (g < b ? 6 : 0);
						break;

					case g:
						hue = (b - r) / delta + 2;
						break;

					case b:
						hue = (r - g) / delta + 4;
						break;
				}

				hue /= 6;
			}

			target.h = hue;
			target.s = saturation;
			target.l = lightness;
			return target;
		};

		_proto.getStyle = function getStyle() {
			return 'rgb(' + (this.r * 255 | 0) + ',' + (this.g * 255 | 0) + ',' + (this.b * 255 | 0) + ')';
		};

		_proto.offsetHSL = function offsetHSL(h, s, l) {
			this.getHSL(_hslA);
			_hslA.h += h;
			_hslA.s += s;
			_hslA.l += l;
			this.setHSL(_hslA.h, _hslA.s, _hslA.l);
			return this;
		};

		_proto.add = function add(color) {
			this.r += color.r;
			this.g += color.g;
			this.b += color.b;
			return this;
		};

		_proto.addColors = function addColors(color1, color2) {
			this.r = color1.r + color2.r;
			this.g = color1.g + color2.g;
			this.b = color1.b + color2.b;
			return this;
		};

		_proto.addScalar = function addScalar(s) {
			this.r += s;
			this.g += s;
			this.b += s;
			return this;
		};

		_proto.sub = function sub(color) {
			this.r = Math.max(0, this.r - color.r);
			this.g = Math.max(0, this.g - color.g);
			this.b = Math.max(0, this.b - color.b);
			return this;
		};

		_proto.multiply = function multiply(color) {
			this.r *= color.r;
			this.g *= color.g;
			this.b *= color.b;
			return this;
		};

		_proto.multiplyScalar = function multiplyScalar(s) {
			this.r *= s;
			this.g *= s;
			this.b *= s;
			return this;
		};

		_proto.lerp = function lerp(color, alpha) {
			this.r += (color.r - this.r) * alpha;
			this.g += (color.g - this.g) * alpha;
			this.b += (color.b - this.b) * alpha;
			return this;
		};

		_proto.lerpHSL = function lerpHSL(color, alpha) {
			this.getHSL(_hslA);
			color.getHSL(_hslB);
			var h = MathUtils.lerp(_hslA.h, _hslB.h, alpha);
			var s = MathUtils.lerp(_hslA.s, _hslB.s, alpha);
			var l = MathUtils.lerp(_hslA.l, _hslB.l, alpha);
			this.setHSL(h, s, l);
			return this;
		};

		_proto.equals = function equals(c) {
			return c.r === this.r && c.g === this.g && c.b === this.b;
		};

		_proto.fromArray = function fromArray(array, offset) {
			if (offset === undefined) offset = 0;
			this.r = array[offset];
			this.g = array[offset + 1];
			this.b = array[offset + 2];
			return this;
		};

		_proto.toArray = function toArray(array, offset) {
			if (array === undefined) array = [];
			if (offset === undefined) offset = 0;
			array[offset] = this.r;
			array[offset + 1] = this.g;
			array[offset + 2] = this.b;
			return array;
		};

		_proto.fromBufferAttribute = function fromBufferAttribute(attribute, index) {
			this.r = attribute.getX(index);
			this.g = attribute.getY(index);
			this.b = attribute.getZ(index);

			if (attribute.normalized === true) {
				// assuming Uint8Array
				this.r /= 255;
				this.g /= 255;
				this.b /= 255;
			}

			return this;
		};

		_proto.toJSON = function toJSON() {
			return this.getHex();
		};

		return Color;
	}();

	Color.NAMES = _colorKeywords;
	Color.prototype.r = 1;
	Color.prototype.g = 1;
	Color.prototype.b = 1;

	/**
	 * Port from https://github.com/mapbox/earcut (v2.2.2)
	 */
	var Earcut = {
		triangulate: function triangulate(data, holeIndices, dim) {
			dim = dim || 2;
			var hasHoles = holeIndices && holeIndices.length;
			var outerLen = hasHoles ? holeIndices[0] * dim : data.length;
			var outerNode = linkedList(data, 0, outerLen, dim, true);
			var triangles = [];
			if (!outerNode || outerNode.next === outerNode.prev) return triangles;
			var minX, minY, maxX, maxY, x, y, invSize;
			if (hasHoles) outerNode = eliminateHoles(data, holeIndices, outerNode, dim); // if the shape is not too simple, we'll use z-order curve hash later; calculate polygon bbox

			if (data.length > 80 * dim) {
				minX = maxX = data[0];
				minY = maxY = data[1];

				for (var i = dim; i < outerLen; i += dim) {
					x = data[i];
					y = data[i + 1];
					if (x < minX) minX = x;
					if (y < minY) minY = y;
					if (x > maxX) maxX = x;
					if (y > maxY) maxY = y;
				} // minX, minY and invSize are later used to transform coords into integers for z-order calculation


				invSize = Math.max(maxX - minX, maxY - minY);
				invSize = invSize !== 0 ? 1 / invSize : 0;
			}

			earcutLinked(outerNode, triangles, dim, minX, minY, invSize);
			return triangles;
		}
	}; // create a circular doubly linked list from polygon points in the specified winding order

	function linkedList(data, start, end, dim, clockwise) {
		var i, last;

		if (clockwise === signedArea(data, start, end, dim) > 0) {
			for (i = start; i < end; i += dim) {
				last = insertNode(i, data[i], data[i + 1], last);
			}
		} else {
			for (i = end - dim; i >= start; i -= dim) {
				last = insertNode(i, data[i], data[i + 1], last);
			}
		}

		if (last && equals(last, last.next)) {
			removeNode(last);
			last = last.next;
		}

		return last;
	} // eliminate colinear or duplicate points


	function filterPoints(start, end) {
		if (!start) return start;
		if (!end) end = start;
		var p = start,
				again;

		do {
			again = false;

			if (!p.steiner && (equals(p, p.next) || area(p.prev, p, p.next) === 0)) {
				removeNode(p);
				p = end = p.prev;
				if (p === p.next) break;
				again = true;
			} else {
				p = p.next;
			}
		} while (again || p !== end);

		return end;
	} // main ear slicing loop which triangulates a polygon (given as a linked list)


	function earcutLinked(ear, triangles, dim, minX, minY, invSize, pass) {
		if (!ear) return; // interlink polygon nodes in z-order

		if (!pass && invSize) indexCurve(ear, minX, minY, invSize);
		var stop = ear,
				prev,
				next; // iterate through ears, slicing them one by one

		while (ear.prev !== ear.next) {
			prev = ear.prev;
			next = ear.next;

			if (invSize ? isEarHashed(ear, minX, minY, invSize) : isEar(ear)) {
				// cut off the triangle
				triangles.push(prev.i / dim);
				triangles.push(ear.i / dim);
				triangles.push(next.i / dim);
				removeNode(ear); // skipping the next vertex leads to less sliver triangles

				ear = next.next;
				stop = next.next;
				continue;
			}

			ear = next; // if we looped through the whole remaining polygon and can't find any more ears

			if (ear === stop) {
				// try filtering points and slicing again
				if (!pass) {
					earcutLinked(filterPoints(ear), triangles, dim, minX, minY, invSize, 1); // if this didn't work, try curing all small self-intersections locally
				} else if (pass === 1) {
					ear = cureLocalIntersections(filterPoints(ear), triangles, dim);
					earcutLinked(ear, triangles, dim, minX, minY, invSize, 2); // as a last resort, try splitting the remaining polygon into two
				} else if (pass === 2) {
					splitEarcut(ear, triangles, dim, minX, minY, invSize);
				}

				break;
			}
		}
	} // check whether a polygon node forms a valid ear with adjacent nodes


	function isEar(ear) {
		var a = ear.prev,
				b = ear,
				c = ear.next;
		if (area(a, b, c) >= 0) return false; // reflex, can't be an ear
		// now make sure we don't have other points inside the potential ear

		var p = ear.next.next;

		while (p !== ear.prev) {
			if (pointInTriangle(a.x, a.y, b.x, b.y, c.x, c.y, p.x, p.y) && area(p.prev, p, p.next) >= 0) return false;
			p = p.next;
		}

		return true;
	}

	function isEarHashed(ear, minX, minY, invSize) {
		var a = ear.prev,
				b = ear,
				c = ear.next;
		if (area(a, b, c) >= 0) return false; // reflex, can't be an ear
		// triangle bbox; min & max are calculated like this for speed

		var minTX = a.x < b.x ? a.x < c.x ? a.x : c.x : b.x < c.x ? b.x : c.x,
				minTY = a.y < b.y ? a.y < c.y ? a.y : c.y : b.y < c.y ? b.y : c.y,
				maxTX = a.x > b.x ? a.x > c.x ? a.x : c.x : b.x > c.x ? b.x : c.x,
				maxTY = a.y > b.y ? a.y > c.y ? a.y : c.y : b.y > c.y ? b.y : c.y; // z-order range for the current triangle bbox;

		var minZ = zOrder(minTX, minTY, minX, minY, invSize),
				maxZ = zOrder(maxTX, maxTY, minX, minY, invSize);
		var p = ear.prevZ,
				n = ear.nextZ; // look for points inside the triangle in both directions

		while (p && p.z >= minZ && n && n.z <= maxZ) {
			if (p !== ear.prev && p !== ear.next && pointInTriangle(a.x, a.y, b.x, b.y, c.x, c.y, p.x, p.y) && area(p.prev, p, p.next) >= 0) return false;
			p = p.prevZ;
			if (n !== ear.prev && n !== ear.next && pointInTriangle(a.x, a.y, b.x, b.y, c.x, c.y, n.x, n.y) && area(n.prev, n, n.next) >= 0) return false;
			n = n.nextZ;
		} // look for remaining points in decreasing z-order


		while (p && p.z >= minZ) {
			if (p !== ear.prev && p !== ear.next && pointInTriangle(a.x, a.y, b.x, b.y, c.x, c.y, p.x, p.y) && area(p.prev, p, p.next) >= 0) return false;
			p = p.prevZ;
		} // look for remaining points in increasing z-order


		while (n && n.z <= maxZ) {
			if (n !== ear.prev && n !== ear.next && pointInTriangle(a.x, a.y, b.x, b.y, c.x, c.y, n.x, n.y) && area(n.prev, n, n.next) >= 0) return false;
			n = n.nextZ;
		}

		return true;
	} // go through all polygon nodes and cure small local self-intersections


	function cureLocalIntersections(start, triangles, dim) {
		var p = start;

		do {
			var a = p.prev,
					b = p.next.next;

			if (!equals(a, b) && intersects(a, p, p.next, b) && locallyInside(a, b) && locallyInside(b, a)) {
				triangles.push(a.i / dim);
				triangles.push(p.i / dim);
				triangles.push(b.i / dim); // remove two nodes involved

				removeNode(p);
				removeNode(p.next);
				p = start = b;
			}

			p = p.next;
		} while (p !== start);

		return filterPoints(p);
	} // try splitting polygon into two and triangulate them independently


	function splitEarcut(start, triangles, dim, minX, minY, invSize) {
		// look for a valid diagonal that divides the polygon into two
		var a = start;

		do {
			var b = a.next.next;

			while (b !== a.prev) {
				if (a.i !== b.i && isValidDiagonal(a, b)) {
					// split the polygon in two by the diagonal
					var c = splitPolygon(a, b); // filter colinear points around the cuts

					a = filterPoints(a, a.next);
					c = filterPoints(c, c.next); // run earcut on each half

					earcutLinked(a, triangles, dim, minX, minY, invSize);
					earcutLinked(c, triangles, dim, minX, minY, invSize);
					return;
				}

				b = b.next;
			}

			a = a.next;
		} while (a !== start);
	} // link every hole into the outer loop, producing a single-ring polygon without holes


	function eliminateHoles(data, holeIndices, outerNode, dim) {
		var queue = [];
		var i, len, start, end, list;

		for (i = 0, len = holeIndices.length; i < len; i++) {
			start = holeIndices[i] * dim;
			end = i < len - 1 ? holeIndices[i + 1] * dim : data.length;
			list = linkedList(data, start, end, dim, false);
			if (list === list.next) list.steiner = true;
			queue.push(getLeftmost(list));
		}

		queue.sort(compareX); // process holes from left to right

		for (i = 0; i < queue.length; i++) {
			eliminateHole(queue[i], outerNode);
			outerNode = filterPoints(outerNode, outerNode.next);
		}

		return outerNode;
	}

	function compareX(a, b) {
		return a.x - b.x;
	} // find a bridge between vertices that connects hole with an outer ring and and link it


	function eliminateHole(hole, outerNode) {
		outerNode = findHoleBridge(hole, outerNode);

		if (outerNode) {
			var b = splitPolygon(outerNode, hole); // filter collinear points around the cuts

			filterPoints(outerNode, outerNode.next);
			filterPoints(b, b.next);
		}
	} // David Eberly's algorithm for finding a bridge between hole and outer polygon


	function findHoleBridge(hole, outerNode) {
		var p = outerNode;
		var hx = hole.x;
		var hy = hole.y;
		var qx = -Infinity,
				m; // find a segment intersected by a ray from the hole's leftmost point to the left;
		// segment's endpoint with lesser x will be potential connection point

		do {
			if (hy <= p.y && hy >= p.next.y && p.next.y !== p.y) {
				var x = p.x + (hy - p.y) * (p.next.x - p.x) / (p.next.y - p.y);

				if (x <= hx && x > qx) {
					qx = x;

					if (x === hx) {
						if (hy === p.y) return p;
						if (hy === p.next.y) return p.next;
					}

					m = p.x < p.next.x ? p : p.next;
				}
			}

			p = p.next;
		} while (p !== outerNode);

		if (!m) return null;
		if (hx === qx) return m; // hole touches outer segment; pick leftmost endpoint
		// look for points inside the triangle of hole point, segment intersection and endpoint;
		// if there are no points found, we have a valid connection;
		// otherwise choose the point of the minimum angle with the ray as connection point

		var stop = m,
				mx = m.x,
				my = m.y;
		var tanMin = Infinity,
				tan;
		p = m;

		do {
			if (hx >= p.x && p.x >= mx && hx !== p.x && pointInTriangle(hy < my ? hx : qx, hy, mx, my, hy < my ? qx : hx, hy, p.x, p.y)) {
				tan = Math.abs(hy - p.y) / (hx - p.x); // tangential

				if (locallyInside(p, hole) && (tan < tanMin || tan === tanMin && (p.x > m.x || p.x === m.x && sectorContainsSector(m, p)))) {
					m = p;
					tanMin = tan;
				}
			}

			p = p.next;
		} while (p !== stop);

		return m;
	} // whether sector in vertex m contains sector in vertex p in the same coordinates


	function sectorContainsSector(m, p) {
		return area(m.prev, m, p.prev) < 0 && area(p.next, m, m.next) < 0;
	} // interlink polygon nodes in z-order


	function indexCurve(start, minX, minY, invSize) {
		var p = start;

		do {
			if (p.z === null) p.z = zOrder(p.x, p.y, minX, minY, invSize);
			p.prevZ = p.prev;
			p.nextZ = p.next;
			p = p.next;
		} while (p !== start);

		p.prevZ.nextZ = null;
		p.prevZ = null;
		sortLinked(p);
	} // Simon Tatham's linked list merge sort algorithm
	// http://www.chiark.greenend.org.uk/~sgtatham/algorithms/listsort.html


	function sortLinked(list) {
		var i,
				p,
				q,
				e,
				tail,
				numMerges,
				pSize,
				qSize,
				inSize = 1;

		do {
			p = list;
			list = null;
			tail = null;
			numMerges = 0;

			while (p) {
				numMerges++;
				q = p;
				pSize = 0;

				for (i = 0; i < inSize; i++) {
					pSize++;
					q = q.nextZ;
					if (!q) break;
				}

				qSize = inSize;

				while (pSize > 0 || qSize > 0 && q) {
					if (pSize !== 0 && (qSize === 0 || !q || p.z <= q.z)) {
						e = p;
						p = p.nextZ;
						pSize--;
					} else {
						e = q;
						q = q.nextZ;
						qSize--;
					}

					if (tail) tail.nextZ = e;else list = e;
					e.prevZ = tail;
					tail = e;
				}

				p = q;
			}

			tail.nextZ = null;
			inSize *= 2;
		} while (numMerges > 1);

		return list;
	} // z-order of a point given coords and inverse of the longer side of data bbox


	function zOrder(x, y, minX, minY, invSize) {
		// coords are transformed into non-negative 15-bit integer range
		x = 32767 * (x - minX) * invSize;
		y = 32767 * (y - minY) * invSize;
		x = (x | x << 8) & 0x00FF00FF;
		x = (x | x << 4) & 0x0F0F0F0F;
		x = (x | x << 2) & 0x33333333;
		x = (x | x << 1) & 0x55555555;
		y = (y | y << 8) & 0x00FF00FF;
		y = (y | y << 4) & 0x0F0F0F0F;
		y = (y | y << 2) & 0x33333333;
		y = (y | y << 1) & 0x55555555;
		return x | y << 1;
	} // find the leftmost node of a polygon ring


	function getLeftmost(start) {
		var p = start,
				leftmost = start;

		do {
			if (p.x < leftmost.x || p.x === leftmost.x && p.y < leftmost.y) leftmost = p;
			p = p.next;
		} while (p !== start);

		return leftmost;
	} // check if a point lies within a convex triangle


	function pointInTriangle(ax, ay, bx, by, cx, cy, px, py) {
		return (cx - px) * (ay - py) - (ax - px) * (cy - py) >= 0 && (ax - px) * (by - py) - (bx - px) * (ay - py) >= 0 && (bx - px) * (cy - py) - (cx - px) * (by - py) >= 0;
	} // check if a diagonal between two polygon nodes is valid (lies in polygon interior)


	function isValidDiagonal(a, b) {
		return a.next.i !== b.i && a.prev.i !== b.i && !intersectsPolygon(a, b) && ( // dones't intersect other edges
		locallyInside(a, b) && locallyInside(b, a) && middleInside(a, b) && ( // locally visible
		area(a.prev, a, b.prev) || area(a, b.prev, b)) || // does not create opposite-facing sectors
		equals(a, b) && area(a.prev, a, a.next) > 0 && area(b.prev, b, b.next) > 0); // special zero-length case
	} // signed area of a triangle


	function area(p, q, r) {
		return (q.y - p.y) * (r.x - q.x) - (q.x - p.x) * (r.y - q.y);
	} // check if two points are equal


	function equals(p1, p2) {
		return p1.x === p2.x && p1.y === p2.y;
	} // check if two segments intersect


	function intersects(p1, q1, p2, q2) {
		var o1 = sign(area(p1, q1, p2));
		var o2 = sign(area(p1, q1, q2));
		var o3 = sign(area(p2, q2, p1));
		var o4 = sign(area(p2, q2, q1));
		if (o1 !== o2 && o3 !== o4) return true; // general case

		if (o1 === 0 && onSegment(p1, p2, q1)) return true; // p1, q1 and p2 are collinear and p2 lies on p1q1

		if (o2 === 0 && onSegment(p1, q2, q1)) return true; // p1, q1 and q2 are collinear and q2 lies on p1q1

		if (o3 === 0 && onSegment(p2, p1, q2)) return true; // p2, q2 and p1 are collinear and p1 lies on p2q2

		if (o4 === 0 && onSegment(p2, q1, q2)) return true; // p2, q2 and q1 are collinear and q1 lies on p2q2

		return false;
	} // for collinear points p, q, r, check if point q lies on segment pr


	function onSegment(p, q, r) {
		return q.x <= Math.max(p.x, r.x) && q.x >= Math.min(p.x, r.x) && q.y <= Math.max(p.y, r.y) && q.y >= Math.min(p.y, r.y);
	}

	function sign(num) {
		return num > 0 ? 1 : num < 0 ? -1 : 0;
	} // check if a polygon diagonal intersects any polygon segments


	function intersectsPolygon(a, b) {
		var p = a;

		do {
			if (p.i !== a.i && p.next.i !== a.i && p.i !== b.i && p.next.i !== b.i && intersects(p, p.next, a, b)) return true;
			p = p.next;
		} while (p !== a);

		return false;
	} // check if a polygon diagonal is locally inside the polygon


	function locallyInside(a, b) {
		return area(a.prev, a, a.next) < 0 ? area(a, b, a.next) >= 0 && area(a, a.prev, b) >= 0 : area(a, b, a.prev) < 0 || area(a, a.next, b) < 0;
	} // check if the middle point of a polygon diagonal is inside the polygon


	function middleInside(a, b) {
		var p = a,
				inside = false;
		var px = (a.x + b.x) / 2,
				py = (a.y + b.y) / 2;

		do {
			if (p.y > py !== p.next.y > py && p.next.y !== p.y && px < (p.next.x - p.x) * (py - p.y) / (p.next.y - p.y) + p.x) inside = !inside;
			p = p.next;
		} while (p !== a);

		return inside;
	} // link two polygon vertices with a bridge; if the vertices belong to the same ring, it splits polygon into two;
	// if one belongs to the outer ring and another to a hole, it merges it into a single ring


	function splitPolygon(a, b) {
		var a2 = new Node(a.i, a.x, a.y),
				b2 = new Node(b.i, b.x, b.y),
				an = a.next,
				bp = b.prev;
		a.next = b;
		b.prev = a;
		a2.next = an;
		an.prev = a2;
		b2.next = a2;
		a2.prev = b2;
		bp.next = b2;
		b2.prev = bp;
		return b2;
	} // create a node and optionally link it with previous one (in a circular doubly linked list)


	function insertNode(i, x, y, last) {
		var p = new Node(i, x, y);

		if (!last) {
			p.prev = p;
			p.next = p;
		} else {
			p.next = last.next;
			p.prev = last;
			last.next.prev = p;
			last.next = p;
		}

		return p;
	}

	function removeNode(p) {
		p.next.prev = p.prev;
		p.prev.next = p.next;
		if (p.prevZ) p.prevZ.nextZ = p.nextZ;
		if (p.nextZ) p.nextZ.prevZ = p.prevZ;
	}

	function Node(i, x, y) {
		// vertex index in coordinates array
		this.i = i; // vertex coordinates

		this.x = x;
		this.y = y; // previous and next vertex nodes in a polygon ring

		this.prev = null;
		this.next = null; // z-order curve value

		this.z = null; // previous and next nodes in z-order

		this.prevZ = null;
		this.nextZ = null; // indicates whether this is a steiner point

		this.steiner = false;
	}

	function signedArea(data, start, end, dim) {
		var sum = 0;

		for (var i = start, j = end - dim; i < end; i += dim) {
			sum += (data[j] - data[i]) * (data[i + 1] + data[j + 1]);
			j = i;
		}

		return sum;
	}

	var ShapeUtils = {
		// calculate area of the contour polygon
		area: function area(contour) {
			var n = contour.length;
			var a = 0.0;

			for (var p = n - 1, q = 0; q < n; p = q++) {
				a += contour[p].x * contour[q].y - contour[q].x * contour[p].y;
			}

			return a * 0.5;
		},
		isClockWise: function isClockWise(pts) {
			return ShapeUtils.area(pts) < 0;
		},
		triangulateShape: function triangulateShape(contour, holes) {
			var vertices = []; // flat array of vertices like [ x0,y0, x1,y1, x2,y2, ... ]

			var holeIndices = []; // array of hole indices

			var faces = []; // final array of vertex indices like [ [ a,b,d ], [ b,c,d ] ]

			removeDupEndPts(contour);
			addContour(vertices, contour); //

			var holeIndex = contour.length;
			holes.forEach(removeDupEndPts);

			for (var i = 0; i < holes.length; i++) {
				holeIndices.push(holeIndex);
				holeIndex += holes[i].length;
				addContour(vertices, holes[i]);
			} //


			var triangles = Earcut.triangulate(vertices, holeIndices); //

			for (var _i = 0; _i < triangles.length; _i += 3) {
				faces.push(triangles.slice(_i, _i + 3));
			}

			return faces;
		}
	};

	function removeDupEndPts(points) {
		var l = points.length;

		if (l > 2 && points[l - 1].equals(points[0])) {
			points.pop();
		}
	}

	function addContour(vertices, contour) {
		for (var i = 0; i < contour.length; i++) {
			vertices.push(contour[i].x);
			vertices.push(contour[i].y);
		}
	}

	function ShapePath() {
		this.type = 'ShapePath';
		this.color = new Color();
		this.subPaths = [];
		this.currentPath = null;
	}

	Object.assign(ShapePath.prototype, {
		moveTo: function moveTo(x, y) {
			this.currentPath = new Path();
			this.subPaths.push(this.currentPath);
			this.currentPath.moveTo(x, y);
			return this;
		},
		lineTo: function lineTo(x, y) {
			this.currentPath.lineTo(x, y);
			return this;
		},
		quadraticCurveTo: function quadraticCurveTo(aCPx, aCPy, aX, aY) {
			this.currentPath.quadraticCurveTo(aCPx, aCPy, aX, aY);
			return this;
		},
		bezierCurveTo: function bezierCurveTo(aCP1x, aCP1y, aCP2x, aCP2y, aX, aY) {
			this.currentPath.bezierCurveTo(aCP1x, aCP1y, aCP2x, aCP2y, aX, aY);
			return this;
		},
		splineThru: function splineThru(pts) {
			this.currentPath.splineThru(pts);
			return this;
		},
		toShapes: function toShapes(isCCW, noHoles) {
			function toShapesNoHoles(inSubpaths) {
				var shapes = [];

				for (var i = 0, l = inSubpaths.length; i < l; i++) {
					var _tmpPath = inSubpaths[i];

					var _tmpShape = new Shape();

					_tmpShape.curves = _tmpPath.curves;
					shapes.push(_tmpShape);
				}

				return shapes;
			}

			function isPointInsidePolygon(inPt, inPolygon) {
				var polyLen = inPolygon.length; // inPt on polygon contour => immediate success		or
				// toggling of inside/outside at every single! intersection point of an edge
				//	with the horizontal line through inPt, left of inPt
				//	not counting lowerY endpoints of edges and whole edges on that line

				var inside = false;

				for (var p = polyLen - 1, q = 0; q < polyLen; p = q++) {
					var edgeLowPt = inPolygon[p];
					var edgeHighPt = inPolygon[q];
					var edgeDx = edgeHighPt.x - edgeLowPt.x;
					var edgeDy = edgeHighPt.y - edgeLowPt.y;

					if (Math.abs(edgeDy) > Number.EPSILON) {
						// not parallel
						if (edgeDy < 0) {
							edgeLowPt = inPolygon[q];
							edgeDx = -edgeDx;
							edgeHighPt = inPolygon[p];
							edgeDy = -edgeDy;
						}

						if (inPt.y < edgeLowPt.y || inPt.y > edgeHighPt.y) continue;

						if (inPt.y === edgeLowPt.y) {
							if (inPt.x === edgeLowPt.x) return true; // inPt is on contour ?
							// continue;				// no intersection or edgeLowPt => doesn't count !!!
						} else {
							var perpEdge = edgeDy * (inPt.x - edgeLowPt.x) - edgeDx * (inPt.y - edgeLowPt.y);
							if (perpEdge === 0) return true; // inPt is on contour ?

							if (perpEdge < 0) continue;
							inside = !inside; // true intersection left of inPt
						}
					} else {
						// parallel or collinear
						if (inPt.y !== edgeLowPt.y) continue; // parallel
						// edge lies on the same horizontal line as inPt

						if (edgeHighPt.x <= inPt.x && inPt.x <= edgeLowPt.x || edgeLowPt.x <= inPt.x && inPt.x <= edgeHighPt.x) return true; // inPt: Point on contour !
						// continue;
					}
				}

				return inside;
			}

			var isClockWise = ShapeUtils.isClockWise;
			var subPaths = this.subPaths;
			if (subPaths.length === 0) return [];
			if (noHoles === true) return toShapesNoHoles(subPaths);
			var solid, tmpPath, tmpShape;
			var shapes = [];

			if (subPaths.length === 1) {
				tmpPath = subPaths[0];
				tmpShape = new Shape();
				tmpShape.curves = tmpPath.curves;
				shapes.push(tmpShape);
				return shapes;
			}

			var holesFirst = !isClockWise(subPaths[0].getPoints());
			holesFirst = isCCW ? !holesFirst : holesFirst; // console.log("Holes first", holesFirst);

			var betterShapeHoles = [];
			var newShapes = [];
			var newShapeHoles = [];
			var mainIdx = 0;
			var tmpPoints;
			newShapes[mainIdx] = undefined;
			newShapeHoles[mainIdx] = [];

			for (var i = 0, l = subPaths.length; i < l; i++) {
				tmpPath = subPaths[i];
				tmpPoints = tmpPath.getPoints();
				solid = isClockWise(tmpPoints);
				solid = isCCW ? !solid : solid;

				if (solid) {
					if (!holesFirst && newShapes[mainIdx]) mainIdx++;
					newShapes[mainIdx] = {
						s: new Shape(),
						p: tmpPoints
					};
					newShapes[mainIdx].s.curves = tmpPath.curves;
					if (holesFirst) mainIdx++;
					newShapeHoles[mainIdx] = []; //console.log('cw', i);
				} else {
					newShapeHoles[mainIdx].push({
						h: tmpPath,
						p: tmpPoints[0]
					}); //console.log('ccw', i);
				}
			} // only Holes? -> probably all Shapes with wrong orientation


			if (!newShapes[0]) return toShapesNoHoles(subPaths);

			if (newShapes.length > 1) {
				var ambiguous = false;
				var toChange = [];

				for (var sIdx = 0, sLen = newShapes.length; sIdx < sLen; sIdx++) {
					betterShapeHoles[sIdx] = [];
				}

				for (var _sIdx = 0, _sLen = newShapes.length; _sIdx < _sLen; _sIdx++) {
					var sho = newShapeHoles[_sIdx];

					for (var hIdx = 0; hIdx < sho.length; hIdx++) {
						var ho = sho[hIdx];
						var hole_unassigned = true;

						for (var s2Idx = 0; s2Idx < newShapes.length; s2Idx++) {
							if (isPointInsidePolygon(ho.p, newShapes[s2Idx].p)) {
								if (_sIdx !== s2Idx) toChange.push({
									froms: _sIdx,
									tos: s2Idx,
									hole: hIdx
								});

								if (hole_unassigned) {
									hole_unassigned = false;
									betterShapeHoles[s2Idx].push(ho);
								} else {
									ambiguous = true;
								}
							}
						}

						if (hole_unassigned) {
							betterShapeHoles[_sIdx].push(ho);
						}
					}
				} // console.log("ambiguous: ", ambiguous);


				if (toChange.length > 0) {
					// console.log("to change: ", toChange);
					if (!ambiguous) newShapeHoles = betterShapeHoles;
				}
			}

			var tmpHoles;

			for (var _i = 0, il = newShapes.length; _i < il; _i++) {
				tmpShape = newShapes[_i].s;
				shapes.push(tmpShape);
				tmpHoles = newShapeHoles[_i];

				for (var j = 0, jl = tmpHoles.length; j < jl; j++) {
					tmpShape.holes.push(tmpHoles[j].h);
				}
			} //console.log("shape", shapes);


			return shapes;
		}
	});

	function Font(data) {
		this.type = 'Font';
		this.data = data;
	}

	Object.assign(Font.prototype, {
		isFont: true,
		generateShapes: function generateShapes(text, size) {
			if (size === undefined) size = 100;
			var shapes = [];
			var paths = createPaths(text, size, this.data);

			for (var p = 0, pl = paths.length; p < pl; p++) {
				Array.prototype.push.apply(shapes, paths[p].toShapes());
			}

			return shapes;
		}
	});

	function createPaths(text, size, data) {
		var chars = Array.from ? Array.from(text) : String(text).split(''); // workaround for IE11, see #13988

		var scale = size / data.resolution;
		var line_height = (data.boundingBox.yMax - data.boundingBox.yMin + data.underlineThickness) * scale;
		var paths = [];
		var offsetX = 0,
				offsetY = 0;

		for (var i = 0; i < chars.length; i++) {
			var char = chars[i];

			if (char === '\n') {
				offsetX = 0;
				offsetY -= line_height;
			} else {
				var ret = createPath(char, scale, offsetX, offsetY, data);
				offsetX += ret.offsetX;
				paths.push(ret.path);
			}
		}

		return paths;
	}

	function createPath(char, scale, offsetX, offsetY, data) {
		var glyph = data.glyphs[char] || data.glyphs['?'];

		if (!glyph) {
			console.error('THREE.Font: character "' + char + '" does not exists in font family ' + data.familyName + '.');
			return;
		}

		var path = new ShapePath();
		var x, y, cpx, cpy, cpx1, cpy1, cpx2, cpy2;

		if (glyph.o) {
			var outline = glyph._cachedOutline || (glyph._cachedOutline = glyph.o.split(' '));

			for (var i = 0, l = outline.length; i < l;) {
				var action = outline[i++];

				switch (action) {
					case 'm':
						// moveTo
						x = outline[i++] * scale + offsetX;
						y = outline[i++] * scale + offsetY;
						path.moveTo(x, y);
						break;

					case 'l':
						// lineTo
						x = outline[i++] * scale + offsetX;
						y = outline[i++] * scale + offsetY;
						path.lineTo(x, y);
						break;

					case 'q':
						// quadraticCurveTo
						cpx = outline[i++] * scale + offsetX;
						cpy = outline[i++] * scale + offsetY;
						cpx1 = outline[i++] * scale + offsetX;
						cpy1 = outline[i++] * scale + offsetY;
						path.quadraticCurveTo(cpx1, cpy1, cpx, cpy);
						break;

					case 'b':
						// bezierCurveTo
						cpx = outline[i++] * scale + offsetX;
						cpy = outline[i++] * scale + offsetY;
						cpx1 = outline[i++] * scale + offsetX;
						cpy1 = outline[i++] * scale + offsetY;
						cpx2 = outline[i++] * scale + offsetX;
						cpy2 = outline[i++] * scale + offsetY;
						path.bezierCurveTo(cpx1, cpy1, cpx2, cpy2, cpx, cpy);
						break;
				}
			}
		}

		return {
			offsetX: glyph.ha * scale,
			path: path
		};
	}

	exports.ArcCurve = ArcCurve;
	exports.CatmullRomCurve3 = CatmullRomCurve3;
	exports.CubicBezierCurve = CubicBezierCurve;
	exports.CubicBezierCurve3 = CubicBezierCurve3;
	exports.Curve = Curve;
	exports.CurvePath = CurvePath;
	exports.EllipseCurve = EllipseCurve;
	exports.Font = Font;
	exports.LineCurve = LineCurve;
	exports.LineCurve3 = LineCurve3;
	exports.Path = Path;
	exports.QuadraticBezierCurve = QuadraticBezierCurve;
	exports.QuadraticBezierCurve3 = QuadraticBezierCurve3;
	exports.Shape = Shape;
	exports.ShapePath = ShapePath;
	exports.ShapeUtils = ShapeUtils;
	exports.SplineCurve = SplineCurve;

	Object.defineProperty(exports, '__esModule', { value: true });

})));
