'use strict';

Object.defineProperty(exports, "__esModule", {
  value: true
});

var _createClass = function () { function defineProperties(target, props) { for (var i = 0; i < props.length; i++) { var descriptor = props[i]; descriptor.enumerable = descriptor.enumerable || false; descriptor.configurable = true; if ("value" in descriptor) descriptor.writable = true; Object.defineProperty(target, descriptor.key, descriptor); } } return function (Constructor, protoProps, staticProps) { if (protoProps) defineProperties(Constructor.prototype, protoProps); if (staticProps) defineProperties(Constructor, staticProps); return Constructor; }; }();

function _classCallCheck(instance, Constructor) { if (!(instance instanceof Constructor)) { throw new TypeError("Cannot call a class as a function"); } }

var ComplexArray = function () {
  function ComplexArray(other) {
    var arrayType = arguments.length > 1 && arguments[1] !== undefined ? arguments[1] : Float32Array;

    _classCallCheck(this, ComplexArray);

    if (other instanceof ComplexArray) {
      // Copy constuctor.
      this.ArrayType = other.ArrayType;
      this.real = new this.ArrayType(other.real);
      this.imag = new this.ArrayType(other.imag);
    } else {
      this.ArrayType = arrayType;
      // other can be either an array or a number.
      this.real = new this.ArrayType(other);
      this.imag = new this.ArrayType(this.real.length);
    }

    this.length = this.real.length;
  }

  _createClass(ComplexArray, [{
    key: 'toString',
    value: function toString() {
      var components = [];

      this.forEach(function (value, i) {
        components.push('(' + value.real.toFixed(2) + ', ' + value.imag.toFixed(2) + ')');
      });

      return '[' + components.join(', ') + ']';
    }
  }, {
    key: 'forEach',
    value: function forEach(iterator) {
      var n = this.length;
      // For gc efficiency, re-use a single object in the iterator.
      var value = Object.seal(Object.defineProperties({}, {
        real: { writable: true }, imag: { writable: true }
      }));

      for (var i = 0; i < n; i++) {
        value.real = this.real[i];
        value.imag = this.imag[i];
        iterator(value, i, n);
      }
    }

    // In-place mapper.

  }, {
    key: 'map',
    value: function map(mapper) {
      var _this = this;

      this.forEach(function (value, i, n) {
        mapper(value, i, n);
        _this.real[i] = value.real;
        _this.imag[i] = value.imag;
      });

      return this;
    }
  }, {
    key: 'conjugate',
    value: function conjugate() {
      return new ComplexArray(this).map(function (value) {
        value.imag *= -1;
      });
    }
  }, {
    key: 'magnitude',
    value: function magnitude() {
      var mags = new this.ArrayType(this.length);

      this.forEach(function (value, i) {
        mags[i] = Math.sqrt(value.real * value.real + value.imag * value.imag);
      });

      return mags;
    }
  }]);

  return ComplexArray;
}();

exports.default = ComplexArray;