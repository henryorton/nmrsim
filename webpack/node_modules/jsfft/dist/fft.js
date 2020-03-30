'use strict';

Object.defineProperty(exports, "__esModule", {
  value: true
});
exports.ComplexArray = undefined;

var _createClass = function () { function defineProperties(target, props) { for (var i = 0; i < props.length; i++) { var descriptor = props[i]; descriptor.enumerable = descriptor.enumerable || false; descriptor.configurable = true; if ("value" in descriptor) descriptor.writable = true; Object.defineProperty(target, descriptor.key, descriptor); } } return function (Constructor, protoProps, staticProps) { if (protoProps) defineProperties(Constructor.prototype, protoProps); if (staticProps) defineProperties(Constructor, staticProps); return Constructor; }; }();

exports.FFT = FFT;
exports.InvFFT = InvFFT;
exports.frequencyMap = frequencyMap;

var _complex_array = require('./complex_array');

var _complex_array2 = _interopRequireDefault(_complex_array);

function _interopRequireDefault(obj) { return obj && obj.__esModule ? obj : { default: obj }; }

function _classCallCheck(instance, Constructor) { if (!(instance instanceof Constructor)) { throw new TypeError("Cannot call a class as a function"); } }

function _possibleConstructorReturn(self, call) { if (!self) { throw new ReferenceError("this hasn't been initialised - super() hasn't been called"); } return call && (typeof call === "object" || typeof call === "function") ? call : self; }

function _inherits(subClass, superClass) { if (typeof superClass !== "function" && superClass !== null) { throw new TypeError("Super expression must either be null or a function, not " + typeof superClass); } subClass.prototype = Object.create(superClass && superClass.prototype, { constructor: { value: subClass, enumerable: false, writable: true, configurable: true } }); if (superClass) Object.setPrototypeOf ? Object.setPrototypeOf(subClass, superClass) : subClass.__proto__ = superClass; }

// Math constants and functions we need.
var PI = Math.PI;
var SQRT1_2 = Math.SQRT1_2;

function FFT(input) {
  return ensureComplexArray(input).FFT();
};

function InvFFT(input) {
  return ensureComplexArray(input).InvFFT();
};

function frequencyMap(input, filterer) {
  return ensureComplexArray(input).frequencyMap(filterer);
};

var ComplexArray = exports.ComplexArray = function (_baseComplexArray) {
  _inherits(ComplexArray, _baseComplexArray);

  function ComplexArray() {
    _classCallCheck(this, ComplexArray);

    return _possibleConstructorReturn(this, (ComplexArray.__proto__ || Object.getPrototypeOf(ComplexArray)).apply(this, arguments));
  }

  _createClass(ComplexArray, [{
    key: 'FFT',
    value: function FFT() {
      return fft(this, false);
    }
  }, {
    key: 'InvFFT',
    value: function InvFFT() {
      return fft(this, true);
    }

    // Applies a frequency-space filter to input, and returns the real-space
    // filtered input.
    // filterer accepts freq, i, n and modifies freq.real and freq.imag.

  }, {
    key: 'frequencyMap',
    value: function frequencyMap(filterer) {
      return this.FFT().map(filterer).InvFFT();
    }
  }]);

  return ComplexArray;
}(_complex_array2.default);

function ensureComplexArray(input) {
  return input instanceof ComplexArray && input || new ComplexArray(input);
}

function fft(input, inverse) {
  var n = input.length;

  if (n & n - 1) {
    return FFT_Recursive(input, inverse);
  } else {
    return FFT_2_Iterative(input, inverse);
  }
}

function FFT_Recursive(input, inverse) {
  var n = input.length;

  if (n === 1) {
    return input;
  }

  var output = new ComplexArray(n, input.ArrayType);

  // Use the lowest odd factor, so we are able to use FFT_2_Iterative in the
  // recursive transforms optimally.
  var p = LowestOddFactor(n);
  var m = n / p;
  var normalisation = 1 / Math.sqrt(p);
  var recursive_result = new ComplexArray(m, input.ArrayType);

  // Loops go like O(n Σ p_i), where p_i are the prime factors of n.
  // for a power of a prime, p, this reduces to O(n p log_p n)
  for (var j = 0; j < p; j++) {
    for (var i = 0; i < m; i++) {
      recursive_result.real[i] = input.real[i * p + j];
      recursive_result.imag[i] = input.imag[i * p + j];
    }
    // Don't go deeper unless necessary to save allocs.
    if (m > 1) {
      recursive_result = fft(recursive_result, inverse);
    }

    var del_f_r = Math.cos(2 * PI * j / n);
    var del_f_i = (inverse ? -1 : 1) * Math.sin(2 * PI * j / n);
    var f_r = 1;
    var f_i = 0;

    for (var _i = 0; _i < n; _i++) {
      var _real = recursive_result.real[_i % m];
      var _imag = recursive_result.imag[_i % m];

      output.real[_i] += f_r * _real - f_i * _imag;
      output.imag[_i] += f_r * _imag + f_i * _real;

      var _ref = [f_r * del_f_r - f_i * del_f_i, f_i = f_r * del_f_i + f_i * del_f_r];
      f_r = _ref[0];
      f_i = _ref[1];
    }
  }

  // Copy back to input to match FFT_2_Iterative in-placeness
  // TODO: faster way of making this in-place?
  for (var _i2 = 0; _i2 < n; _i2++) {
    input.real[_i2] = normalisation * output.real[_i2];
    input.imag[_i2] = normalisation * output.imag[_i2];
  }

  return input;
}

function FFT_2_Iterative(input, inverse) {
  var n = input.length;

  var output = BitReverseComplexArray(input);
  var output_r = output.real;
  var output_i = output.imag;
  // Loops go like O(n log n):
  //   width ~ log n; i,j ~ n
  var width = 1;
  while (width < n) {
    var del_f_r = Math.cos(PI / width);
    var del_f_i = (inverse ? -1 : 1) * Math.sin(PI / width);
    for (var i = 0; i < n / (2 * width); i++) {
      var f_r = 1;
      var f_i = 0;
      for (var j = 0; j < width; j++) {
        var l_index = 2 * i * width + j;
        var r_index = l_index + width;

        var left_r = output_r[l_index];
        var left_i = output_i[l_index];
        var right_r = f_r * output_r[r_index] - f_i * output_i[r_index];
        var right_i = f_i * output_r[r_index] + f_r * output_i[r_index];

        output_r[l_index] = SQRT1_2 * (left_r + right_r);
        output_i[l_index] = SQRT1_2 * (left_i + right_i);
        output_r[r_index] = SQRT1_2 * (left_r - right_r);
        output_i[r_index] = SQRT1_2 * (left_i - right_i);

        var _ref2 = [f_r * del_f_r - f_i * del_f_i, f_r * del_f_i + f_i * del_f_r];
        f_r = _ref2[0];
        f_i = _ref2[1];
      }
    }
    width <<= 1;
  }

  return output;
}

function BitReverseIndex(index, n) {
  var bitreversed_index = 0;

  while (n > 1) {
    bitreversed_index <<= 1;
    bitreversed_index += index & 1;
    index >>= 1;
    n >>= 1;
  }
  return bitreversed_index;
}

function BitReverseComplexArray(array) {
  var n = array.length;
  var flips = new Set();

  for (var i = 0; i < n; i++) {
    var r_i = BitReverseIndex(i, n);

    if (flips.has(i)) continue;

    var _ref3 = [array.real[r_i], array.real[i]];
    array.real[i] = _ref3[0];
    array.real[r_i] = _ref3[1];
    var _ref4 = [array.imag[r_i], array.imag[i]];
    array.imag[i] = _ref4[0];
    array.imag[r_i] = _ref4[1];


    flips.add(r_i);
  }

  return array;
}

function LowestOddFactor(n) {
  var sqrt_n = Math.sqrt(n);
  var factor = 3;

  while (factor <= sqrt_n) {
    if (n % factor === 0) return factor;
    factor += 2;
  }
  return n;
}