'use strict';

var _fft = require('./fft');

var fft = _interopRequireWildcard(_fft);

function _interopRequireWildcard(obj) { if (obj && obj.__esModule) { return obj; } else { var newObj = {}; if (obj != null) { for (var key in obj) { if (Object.prototype.hasOwnProperty.call(obj, key)) newObj[key] = obj[key]; } } newObj.default = obj; return newObj; } }

fft.FFTImageDataRGBA = function (data, nx, ny) {
  var rgb = splitRGB(data);

  return mergeRGB(FFT2D(new fft.ComplexArray(rgb[0], Float32Array), nx, ny), FFT2D(new fft.ComplexArray(rgb[1], Float32Array), nx, ny), FFT2D(new fft.ComplexArray(rgb[2], Float32Array), nx, ny));
};

function splitRGB(data) {
  var n = data.length / 4;
  var r = new Uint8ClampedArray(n);
  var g = new Uint8ClampedArray(n);
  var b = new Uint8ClampedArray(n);

  for (var i = 0; i < n; i++) {
    r[i] = data[4 * i];
    g[i] = data[4 * i + 1];
    b[i] = data[4 * i + 2];
  }

  return [r, g, b];
}

function mergeRGB(r, g, b) {
  var n = r.length;
  var output = new fft.ComplexArray(n * 4);

  for (var i = 0; i < n; i++) {
    output.real[4 * i] = r.real[i];
    output.imag[4 * i] = r.imag[i];
    output.real[4 * i + 1] = g.real[i];
    output.imag[4 * i + 1] = g.imag[i];
    output.real[4 * i + 2] = b.real[i];
    output.imag[4 * i + 2] = b.imag[i];
  }

  return output;
}

function FFT2D(input, nx, ny, inverse) {
  var transform = inverse ? 'InvFFT' : 'FFT';
  var output = new fft.ComplexArray(input.length, input.ArrayType);
  var row = new fft.ComplexArray(nx, input.ArrayType);
  var col = new fft.ComplexArray(ny, input.ArrayType);

  var _loop = function _loop(j) {
    row.map(function (v, i) {
      v.real = input.real[i + j * nx];
      v.imag = input.imag[i + j * nx];
    });
    row[transform]().forEach(function (v, i) {
      output.real[i + j * nx] = v.real;
      output.imag[i + j * nx] = v.imag;
    });
  };

  for (var j = 0; j < ny; j++) {
    _loop(j);
  }

  var _loop2 = function _loop2(i) {
    col.map(function (v, j) {
      v.real = output.real[i + j * nx];
      v.imag = output.imag[i + j * nx];
    });
    col[transform]().forEach(function (v, j) {
      output.real[i + j * nx] = v.real;
      output.imag[i + j * nx] = v.imag;
    });
  };

  for (var i = 0; i < nx; i++) {
    _loop2(i);
  }

  return output;
}