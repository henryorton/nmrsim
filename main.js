

class Spin {
  constructor(amplitude=1.0, offset=2*math.pi, phase=0.0, r1=1.0, r2=1.0) {
    this.amplitude = amplitude;
    this.offset = offset;
    this.phase = phase;
    this.r1 = r1;
    this.r2 = r2;
  }
  exponent() {
    return math.complex(-this.r2, this.offset - 2*math.PI*carrierFrequency);
  }
}

class AcquisitionParameters {
  constructor() {
    this.lastSet = "numberPoints";
    this.lastLastSet = "acquisitionTime";
    this.numberPoints = 128;
    this.acquisitionTime = 1.0;
    this.variables = ["numberPoints","acquisitionTime","spectralWidth","dwellTime"];

    this.variables.forEach( v => {
      $(`#${v}`).on("focusout", event => {
        this.setVariable(v, $(event.target).val())
      })
    })
  }
  get dwellTime () {
    return this.acquisitionTime / this.numberPoints;
  }
  get spectralWidth () {
    return 1.0 / this.dwellTime;
  }
  validate (variable, revert) {
    if (this.lastSet==variable) {
      if (this.lastLastSet==variable) {
        this.lastLastSet = revert;
      }
      this.lastSet = this.lastLastSet;
    }
  }
  storeCurrent () {
    this.variables.forEach( v => {
      if ( $(`#${v}`).is(":focus") ) {
        this.setVariable(v, $(`#${v}`).val())
      }
    })
  }
  populate () {
    $("#numberPoints").val(this.numberPoints);
    $("#acquisitionTime").val(this.acquisitionTime);
    $("#spectralWidth").val(this.spectralWidth);
    $("#dwellTime").val(this.dwellTime); 
  }
  setVariable (variable, value) {
    switch (variable) {
      case "numberPoints":
      this.setNumberPoints(value);
      break;
      case "acquisitionTime":
      this.setAcquisitionTime(value);
      break;
      case "spectralWidth":
      this.setSpectralWidth(value);
      break;
      case "dwellTime":
      this.setDwellTime(value);
      break;
    }
  }
  setNumberPoints (value) {
    this.validate("numberPoints", "acquisitionTime");

    if (this.lastSet!="acquisitionTime") {
      this.acquisitionTime = this.dwellTime * value;
    }
    this.numberPoints = math.round(value);

    this.lastLastSet = this.lastSet;
    this.lastSet = "numberPoints";
    this.populate();
  }
  setAcquisitionTime (value) {
    this.validate("acquisitionTime", "numberPoints");

    if (this.lastSet!="numberPoints") {
      this.numberPoints = math.round(this.acquisitionTime / this.dwellTime);
    }
    this.acquisitionTime = value;

    this.lastLastSet = this.lastSet;
    this.lastSet = "acquisitionTime"; 
    this.populate();
  }
  setSpectralWidth (value) {
    this.validate("spectralWidth", "acquisitionTime");
    if (this.lastSet=="numberPoints") {
      this.acquisitionTime = this.numberPoints / value;
    }
    else {
      this.numberPoints = math.round(this.acquisitionTime * value);
    }

    this.lastLastSet = this.lastSet;
    this.lastSet = "spectralWidth";
    this.populate();
  }
  setDwellTime (value) {
    this.setSpectralWidth(1.0/value);
  }
}

var spins = new Array(2).fill().map(() => new Spin())
var aPars = new AcquisitionParameters();
var numberPoints;
var acquisitionTime;
var spectralWidth;
var carrierFrequency;
var timeAxis;
var freqAxis;
var timeDomain;
var freqDomain;
var fftdata;

function updateNumberPoint () {
  parseFloat($("#numberPoints").val());
}

function populateVariables () {
  numberPoints = parseInt($("#numberPoints").val());
  acquisitionTime = parseFloat($("#acquisitionTime").val());
  spectralWidth = parseFloat($("#spectralWidth").val());
  carrierFrequency = parseFloat($("#carrierFrequency").val());
  spins.forEach( (spin, i) => {
    spin.amplitude = parseFloat($(`.spins .spin.${i} .amplitude`).val())
    spin.offset = 2.0 * math.PI * parseFloat($(`.spins .spin.${i} .offset`).val())
    spin.r2 = 1.0 / parseFloat($(`.spins .spin.${i} .t2`).val())
    spin.phase = 2.0 * math.PI * parseFloat($(`.spins .spin.${i} .phase`).val())
  })
}

function initAxes () {
  timeAxis = new Array(numberPoints);
  freqAxis = new Array(numberPoints);
  for (let i=0; i<numberPoints; i++) {
    timeAxis[i] = acquisitionTime * (i / numberPoints)
    freqAxis[i] = (0.5 - (i / numberPoints)) * spectralWidth + carrierFrequency
  }
}

function calculateSignal () {
  timeDomain = new Array(numberPoints).fill(math.complex(0.0));
  spins.forEach(
    spin => {
      timeDomain = math.chain(timeAxis)
        .multiply(spin.exponent())
        .exp()
        .multiply(spin.amplitude)
        .add(timeDomain)
        .done()
    }
  )
}

function processSignal () {
  fftdata = new ComplexArray(numberPoints)
  fftdata.real = math.re(timeDomain);
  fftdata.imag = math.im(timeDomain);
  fftdata.FFT();
  fftdata.fftShift();
  freqDomain = new Array(numberPoints)
  fftdata.forEach(
    (value, i, n) => freqDomain[i] = math.complex(value.real, value.imag)
  )
}

var animationDuration = 800;
var tdPlot = initPlot("#plot-time-domain-svg");
var fqPlot = initPlot("#plot-freq-domain-svg");



function initPlot (canvas_id, axis) {

  var canvas = d3.select(canvas_id);
  var svg = canvas.append("g").attr("class", "frame")
  var margin = {top: 20, right: 20, bottom: 40, left: 40};
  var width = canvas.attr("width") - margin.left - margin.right;
  var height = canvas.attr("height") - margin.top - margin.bottom;

  svg.attr("transform", `translate(${margin.left}, ${margin.top})`);

  var xScale = d3.scaleLinear()
    .range([0, width]);

  var yScale = d3.scaleLinear()
    .range([height, 0]);

  svg.append("g")
    .attr("class", "x_axis")
    .attr("transform", "translate(0," + height + ")")

  svg.append("g")
    .attr("class", "y_axis")

  svg.append("path")
    .attr("class", "line_real")

  svg.append("path")
    .attr("class", "line_imag")

  var plotData = {
    canvas : canvas,
    svg: svg,
    margin: margin,
    width: width,
    height: height,
    xScale: xScale,
    yScale: yScale,
    real: svg.select(".line_real"),
    imag: svg.select(".line_imag"),
    xaxis: svg.select(".x_axis"),
    yaxis: svg.select(".y_axis"),
  }

  return plotData
}

function updateTimeDomainPlot (plot) {

  let range = math.max([math.abs(math.re(timeDomain)), math.abs(math.im(timeDomain))])
  plot.yScale.domain([-range, range])
  plot.yaxis.transition()
    .duration(animationDuration)
    .call(d3.axisLeft(plot.yScale))

  plot.xScale.domain([0, acquisitionTime])
  plot.xaxis.transition()
    .duration(animationDuration)
    .call(d3.axisBottom(plot.xScale))

  let line = d3.line()
    .x( (d, i) => plot.xScale(timeAxis[i]) )
    .y(  d     => plot.yScale(d) );

  plot.real.datum(math.re(timeDomain))
    .transition()
    .duration(animationDuration)
    .attr("d", line)

  plot.imag.datum(math.im(timeDomain))
    .transition()
    .duration(animationDuration)
    .attr("d", line)

}


function updateFreqDomainPlot (plot) {

  plot.yScale.domain(d3.extent(math.re(freqDomain).concat(math.im(freqDomain))))
  plot.yaxis.transition()
    .duration(animationDuration)
    .call(d3.axisLeft(plot.yScale))

  plot.xScale.domain([0.5*spectralWidth + carrierFrequency, -0.5*spectralWidth + carrierFrequency])
  plot.xaxis.transition()
    .duration(animationDuration)
    .call(d3.axisBottom(plot.xScale))

  let line = d3.line()
    .x( (d, i) => plot.xScale(freqAxis[i]) )
    .y(  d     => plot.yScale(d) );

  plot.real.datum(math.re(freqDomain))
    .transition()
    .duration(animationDuration)
    .attr("d", line)

  plot.imag.datum(math.im(freqDomain))
    .transition()
    .duration(animationDuration)
    .attr("d", line)

}

function update () {
  aPars.storeCurrent();
  populateVariables();
  initAxes();
  calculateSignal();
  processSignal();
  updateTimeDomainPlot(tdPlot);
  updateFreqDomainPlot(fqPlot);
}


(function($) {
  $.fn.inputFilter = function(inputFilter) {
    return this.on("input keydown keyup mousedown mouseup select contextmenu drop", function() {
      if (inputFilter(this.value)) {
        this.oldValue = this.value;
        this.oldSelectionStart = this.selectionStart;
        this.oldSelectionEnd = this.selectionEnd;
      } else if (this.hasOwnProperty("oldValue")) {
        this.value = this.oldValue;
        this.setSelectionRange(this.oldSelectionStart, this.oldSelectionEnd);
      } else {
        this.value = "";
      }
    });
  };
}(jQuery));


function initInputFilters () {
  $("#numberPoints").inputFilter( value => /^\d*$/.test(value) );
  $("#acquisitionTime").inputFilter( value => /^\d*[.]?\d*$/.test(value) );
  $("#spectralWidth").inputFilter( value => /^\d*[.]?\d*$/.test(value) );
  $("#dwellTime").inputFilter( value => /^\d*[.]?\d*$/.test(value) );
  $("#carrierFrequency").inputFilter( value => /^-?\d*[.]?\d*$/.test(value) );

  $(".spins .amplitude").each( (i, elem) => {
    $(elem).inputFilter( value => /^\d*[.]?\d*$/.test(value) )
  })
  $(".spins .offset").each( (i, elem) => {
    $(elem).inputFilter( value => /^-?\d*[.]?\d*$/.test(value) )
  })
  $(".spins .t2").each( (i, elem) => {
    $(elem).inputFilter( value => /^\d*[.]?\d*$/.test(value) )
  })
  $(".spins .phase").each( (i, elem) => {
    $(elem).inputFilter( value => /^-?\d*[.]?\d*$/.test(value) )
  })
}


$(document).ready( () => {
  $(document).on('keypress',function(e) {
    if(e.which == 13) {
      update();
    }
  });

  initInputFilters();
  update();

 

});


































