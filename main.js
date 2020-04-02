

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
  evolution(time) {
    let mag = math.chain(time)
      .multiply(this.exponent())
      .add(math.complex(0.0, this.phase))
      .exp()
      .multiply(this.amplitude)
      .done()
    return mag
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
var experimentalNoise;
var quadratureDetection;

var zeroFilling;

var timeAxis;
var freqAxis;
var timeDomain;
var freqDomain;
var fftdata;

function updateNumberPoint () {
  parseFloat($("#numberPoints").val());
}

function populateVariables () {
  zeroFilling = parseInt($("#zeroFilling").val());
  numberPoints = parseInt($("#numberPoints").val());
  acquisitionTime = parseFloat($("#acquisitionTime").val());
  spectralWidth = parseFloat($("#spectralWidth").val());
  carrierFrequency = parseFloat($("#carrierFrequency").val());
  experimentalNoise = parseFloat($("#experimentalNoise").val());
  quadratureDetection = $("#quadratureDetection").is($(":checked"));

  spins.forEach( (spin, i) => {
    spin.amplitude = parseFloat($(`.spins .spin${i} .amplitude`).val())
    spin.offset = 2.0 * math.PI * parseFloat($(`.spins .spin${i} .offset`).val())
    spin.r2 = 1.0 / parseFloat($(`.spins .spin${i} .t2`).val())
    spin.phase = (math.PI / 180.0) * parseFloat($(`.spins .spin${i} .phase`).val())
  })
}

function initAxes () {
  timeAxis = new Array(numberPoints);
  freqAxis = new Array(numberPoints);
  for (let i=0; i<numberPoints; i++) {
    timeAxis[i] = acquisitionTime * (i / numberPoints)
    freqAxis[i] = (0.5 - (i / numberPoints)) * spectralWidth + carrierFrequency
  }
  noiseAxis = math.multiply(math.complex(0.5,0.5), 
    math.random([numberPoints], -experimentalNoise, experimentalNoise)
  )

    
}

function calculateSignal () {
  timeDomain = new Array(numberPoints).fill(math.complex(0.0));
  spins.forEach(
    spin => {
      timeDomain = math.chain(timeAxis)
        .multiply(spin.exponent())
        .add(math.complex(0.0, spin.phase))
        .exp()
        .multiply(spin.amplitude)
        .add(timeDomain)
        .add(noiseAxis)
        .done()
    }
  )
  if (!quadratureDetection) {
    timeDomain.forEach( (elem, i) => {
      timeDomain[i].im = 0;
    })
  }
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

var animationDuration = 200;
var tdPlot = initPlot("#plot-time-domain-svg");
var fqPlot = initPlot("#plot-freq-domain-svg");
var vecPlot = initVectorPlot();
initMouseOverEffects(tdPlot);



function initPlot (canvas_id, axis) {

  var canvas = d3.select(canvas_id);
  var svg = canvas.append("g").attr("class", "frame")
  var margin = {top: 20, right: 20, bottom: 30, left: 30};
  var width = canvas.attr("width") - margin.left - margin.right;
  var height = canvas.attr("height") - margin.top - margin.bottom;

  svg.attr("transform", `translate(${margin.left}, ${margin.top})`);

  var xScale = d3.scaleLinear()
    .range([0, width]);

  var yScale = d3.scaleLinear()
    .range([height, 0]);

  svg.append("g")
    .attr("class", "x_axis")
    .attr("transform", `translate(0,${height})`)

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
    range: null
  }

  return plotData
}





function initMouseOverEffects (plot) {
  var mouseG = plot.svg.append("g")
      .attr("class", "mouse-over-effects");

  mouseG.append("path") // this is the black vertical line to follow mouse
    .attr("class", "mouse-line")
    .style("stroke", "black")
    .style("stroke-width", "1px")
    .attr("visibility", "hidden")

  mouseG.append('svg:rect') // append a rect to catch mouse movements on canvas
    .attr('width', plot.width) // can't catch mouse events on a g element
    .attr('height', plot.height)
    .attr('fill', 'none')
    .attr('pointer-events', 'all')
    .on('mouseout', () => { // on mouse out hide line, circles and text
      d3.select(".mouse-line")
        .attr("visibility", "hidden");
      d3.selectAll("#plot-spin-vector-svg .vector")
        .attr("visibility", "hidden");
    })
    .on('mouseover', () => { // on mouse in show line, circles and text
      d3.select(".mouse-line")
        .attr("visibility", "visible");
      d3.selectAll("#plot-spin-vector-svg .vector.sum")
        .attr("visibility", "visible");
      spins.forEach( (spin, i) => {
        if (spin.amplitude) {
          d3.select(`#plot-spin-vector-svg .vector.spin${i}`)
            .attr("visibility", "visible")
        }
      })
    })
    .on('mousemove', () => { // mouse moving over canvas
      let mouse = d3.mouse(d3.event.currentTarget);
      let point = Math.round( (mouse[0]/plot.width) * numberPoints)
      let time = timeAxis[point];
      let loc = plot.xScale(time);
      let mag = timeDomain[point];
      if ( typeof loc == 'number') {
        d3.select(".mouse-line")
          .attr("d", () => `M${loc},${plot.height} ${loc},0`);
        d3.select("#plot-spin-vector-svg .vector.sum")
          .select("line")
          .attr("x2", (mag.re / plot.range) * vecPlot.radius)
          .attr("y2", -(mag.im / plot.range) * vecPlot.radius)
        spins.forEach( (spin, i) => {
          let mag = spin.evolution(time);
          d3.select(`#plot-spin-vector-svg .vector.spin${i}`)
            .select("line")
            .attr("x2", (mag.re / plot.range) * vecPlot.radius)
            .attr("y2", -(mag.im / plot.range) * vecPlot.radius);
          d3.select(`#plot-spin-vector-svg .vector.spin${i} .pointer`)
            .attr("transform", `rotate(${(-180.0/math.PI)*math.arg(mag)})`)
        })
      }
    })
}

function appendArrowheadMarker (svg, class_name) {
  let id = `${class_name}_arrowhead`
  svg.append("svg:defs").append("svg:marker")
    .attr("id", id)
    .attr("class", `${class_name} arrowhead`)
    .attr("refX", 11)
    .attr("refY", 6)
    .attr("markerWidth", 30)
    .attr("markerHeight", 30)
    .attr("orient", "auto")
    .append("path")
    .attr("d", "M 0 0 12 6 0 12 3 6")
  return `url(#${id})`;
}

function test (v) {
  console.log(v);
}

function initVectorPlot () {

  var canvas = d3.select("#plot-spin-vector-svg");
  var svg = canvas.append("g").attr("class", "frame")
  var margin = {top: 40, right: 40, bottom: 40, left: 40};
  var width = canvas.attr("width") - margin.left - margin.right;
  var height = canvas.attr("height") - margin.top - margin.bottom;
  var radius = 0.5*height;
  
  svg.attr("transform", `translate(${margin.left+0.5*width}, ${margin.top+0.5*height})`);

  svg.append("circle")
    .attr("class", "circle")
    .attr("r", radius)
    .attr("fill", "none");
      
  svg.append("g")
    .attr("class", "vector sum")
    .attr("visibility", "hidden")
    .append("line")
    .attr("marker-end", appendArrowheadMarker(svg, "sum"));

  svg.append("g")
    .attr("class", "line_real")
    .append("line")
    .attr("x2", radius+25)
    .attr("marker-end", appendArrowheadMarker(svg, "line_real"));

  svg.append("g")
    .attr("class", "line_imag")
    .append("line")
    .attr("y2", -radius-25)
    .attr("marker-end", appendArrowheadMarker(svg, "line_imag"));

  spins.forEach( (spin, i) => {
    let spinG = svg.append("g")
      .attr("class", `vector spin spin${i}`)
      .attr("visibility", "hidden")
    spinG.append("line")
      .attr("marker-end", appendArrowheadMarker(svg, "spin"))

    let pointerG = spinG.append("g")
      .attr("class", "pointer")
    pointerG.append("path")
      .attr("d", `M ${radius} 0 ${15+radius} 10 ${15+radius} -10`)
    pointerG.append("g")
      .attr("class", "text")
      .append("text")
      .text(`${'I'.repeat(i+1)}`)
      .style("writing-mode", "tb")
      .attr("x", `${radius+10}`)
      .attr("y", "4.5")
      .attr("dominant-baseline", "center")
      .attr("text-anchor", "middle")
      .attr("font-size", "10")
      .attr("fill", "white")
      .attr("stroke", "none")

  })

  return {
    svg: svg,
    radius: radius
  }
}









function updateTimeDomainPlot (plot) {

  let range = math.max([math.abs(math.re(timeDomain)), math.abs(math.im(timeDomain))])
  plot.range = range;
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


































