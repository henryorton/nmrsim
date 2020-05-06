console.log(window.location.href);

// Global variables
const ArrayType = Float64Array;
const MAX_POINTS = (2**12) // 4k data
const MAX_ZEROFILL = (2**12) // 4k zerofill
const MAX_BYTES = (MAX_POINTS + MAX_ZEROFILL) * ArrayType.BYTES_PER_ELEMENT * 2;
const NUMBER_SPINS = 2;
var ANI_DUR = 200;
const PROC_PLOT_POINTS = 200;

const signalDataBuffer = new ArrayBuffer(MAX_BYTES);
const procDataBuffer = new ArrayBuffer(MAX_BYTES);


$(document).ready( () => {
  initAll();
  $(document).on("keypress", e => {
    if (e.which==13 && !$(e.target).hasClass("parameter")) {
      update();
    }
  })
});


function initSignalArray(acquisitionPoints) {
  let n = acquisitionPoints;
  if (n > MAX_POINTS) {
    throw "Number of points too large"
    return false;
  }
  sig = {
    real: new ArrayType(signalDataBuffer, 0, n),
    imag: new ArrayType(signalDataBuffer, MAX_BYTES/2, n),
  };
  sig.length = sig.real.length;
};


function initTimeArray(zeroFillingPoints) {
  let n = sig.length;
  let z = zeroFillingPoints;
  if (z > MAX_ZEROFILL) {
    throw "Number of zero-filling points too large"
    return false;
  }
  time = {
    real: new ArrayType(signalDataBuffer, 0, n+z),
    imag: new ArrayType(signalDataBuffer, MAX_BYTES/2, n+z),
    length: n+z,
  }
  time.real.fill(0.0, n, n+z);
  time.imag.fill(0.0, n, n+z);
}


function initProcessArray() {
  let nz = time.length;
  let procR = new ArrayType(procDataBuffer, 0, nz);
  let procI = new ArrayType(procDataBuffer, MAX_BYTES/2, nz);
  proc = new ComplexArray(0, ArrayType).from_ri(procR, procI);
}


function initSpins(number) {
  let spins = []
  for (let i=0; i<number; i++) {
    spins[i] = new Spin(`.spin${i}`)
  }
  return spins
}


// Regex input filter for text
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


// Parameter class to handle input filters, ranges and types
// Parameter type can be spin, acquisition, processing
class Parameter {
  constructor(input_id, parameter_type, input_type=null, factor=1.0, range=null, bounds=null) {
    this.val = null;
    this.bounds = bounds;
    this.range = range;
    this.factor = factor;
    this.parameter_type = parameter_type;
    this.input_type = input_type;
    this.id = input_id;
    this.elem = $(input_id)
    this.elem.addClass("parameter")
    this.on_store = () => {};

    if (parameter_type=="spin") {
      let s = input_id.split(" ")
      this.check = $(`#${s[0].slice(1)}-${s[1].slice(1)}Check`)
    }
    else {
      this.check = $(`${input_id}Check`)
    }

    this.check.on("change", e => {
      slider.change_parameter(this);
    })

    if (input_type) {
      switch (input_type) {
        case "float":
          this.regex = /^-?\d*[.]?\d*$/;
          this.parser = parseFloat;
          break;
        case "ufloat":
          this.regex = /^\d*[.]?\d*$/;
          this.parser = parseFloat;
          break;
        case "int":
          this.regex = /^-?\d*$/;
          this.parser = parseInt;
          break;
        case "uint":
          this.regex = /^\d*$/;
          this.parser = parseInt;
          break;
        default:
          throw "Regex not recognised"
      }
      this.elem.inputFilter( value => this.regex.test(value) )
    }

    this.store();

    this.elem.on("keypress", e => {
      if (e.which==13) {
        if (this.store()) {
          update(this.parameter_type);
        }
      }
    });

    this.elem.on("focusout", e => this.store());

  }

  store(value=this.elem.val(), onStore=true) {
    let val = this.parser(value);
    if (isNaN(val)) {
      alert("Invalid input");
      this.set_field();
      return false;
    }
    if (this.bounds) {
      if (!(val >= this.bounds[0] && val <= this.bounds[1])) {
        alert(`Input out of range: [${this.bounds}]`);
        this.set_field();
        return false;
      }
    }
    this.val = val * this.factor;
    if (slider.par==this) {
      slider.parameter_value_change()
    }
    if (onStore) {
      this.on_store();
    }
    return true
  }

  set_field(value=this.val) {
    this.elem.val(value / this.factor);
  }

  set_value(value) {
    if (this.store(value)) {
      this.set_field();
      return true;
    }
    else {
      return false;
    }
  }
}


class Spin {
  constructor(spin_class, amplitude=0.0, offset=0.0, phase=0.0, r2=0.0) {
    this.amplitude = new Parameter(`${spin_class} .amplitude`, "spin", "float", 1.0, [-2,2]);
    this.offset = new Parameter(`${spin_class} .offset`, "spin", "float", 2*PI, [-50,50]);
    this.phase = new Parameter(`${spin_class} .phase`, "spin", "float", PI/180, [0,360]);
    this.t2 = new Parameter(`${spin_class} .t2`, "spin", "ufloat", 1.0, [0,1]);
  }

  evolution(time) {
    let pf = this.amplitude.val * Math.exp(-time / this.t2.val);
    let real = pf * Math.cos((this.offset.val * time) + this.phase.val);
    let imag = pf * Math.sin((this.offset.val * time) + this.phase.val);
    return {real:real, imag:imag}
  }
}


class Slider {
  constructor (slider_id, min_id, max_id, points=200) {
    this.slider = $(slider_id);
    this.min = $(min_id)
    this.max = $(max_id)
    this.points = points;
    this.par = null;

    this.slider.on("input", e => {
      this.slider_change()
    })

  }

  get minv () {
    return this.par.parser(this.min.val())
  }

  get maxv () {
    return this.par.parser(this.max.val())
  }

  pos_from_par (value=this.par.val / this.par.factor) {
    return (value - this.minv) * (this.points / (this.maxv - this.minv))
  }

  par_from_pos (value=this.slider.val()) {
    let val = this.par.parser(value);
    return ((val / this.points) * (this.maxv - this.minv)) + this.minv
  }

  change_parameter (parameter) {
    this.par = parameter;
    let val = this.par.val
    let min = this.par.range[0]
    let max = this.par.range[1]
    if (val < min) {
      this.min.val(val)
      this.par.range[0] = val;
    }
    else {
      this.min.val(min)
    }
    if (val > max) {
      this.max.val(val)
      this.par.range[1] = val;
    }
    else {
      this.max.val(max)
    }
    this.slider.val(this.pos_from_par());

    this.min.off();
    this.max.off();
    this.min.inputFilter(value => this.par.regex.test(value))
    this.max.inputFilter(value => this.par.regex.test(value))
    this.min.on("keypress focusout", e => {
      if (e.which == 13 | e.type == "focusout") {
        this.min_change();
      }
    })
    this.max.on("keypress focusout", e => {
      if (e.which == 13 | e.type == "focusout") {
        this.max_change();
      }
    })


  }

  slider_change () {
    this.par.store(this.par_from_pos());
    this.par.set_field();
    let lastAniDur = ANI_DUR;
    ANI_DUR = 0;
    update(this.par.parameter_type)
    ANI_DUR = lastAniDur;
  }

  parameter_value_change () {
    let val = this.par.val / this.par.factor
    if (val < this.minv) {
      this.min.val(val)
      this.par.range[0] = val;
    }
    if (val > this.maxv) {
      this.max.val(val)
      this.par.range[1] = val;
    }
    this.slider.val(this.pos_from_par());
  }

  min_change () {
    var val = this.par.val / this.par.factor;
    var min = this.minv;
    if (isNaN(min)) {
      alert("Invalid input");
      this.change_parameter(this.par);
      return;
    }
    if (this.par.bounds) {
      if (!(min >= this.par.bounds[0] && min <= this.par.bounds[1])) {
        alert(`Input out of range: [${this.par.bounds}]`);
        this.change_parameter(this.par);
        return;
      }
    }
    if (this.minv > this.maxv) {
      this.max.val(this.min.val());
    }
    if (val < this.minv) {
      this.par.store(this.minv)
      this.par.set_field();
      var pos = this.pos_from_par(this.minv);
    }
    else {
      var pos = this.pos_from_par()
    }
    this.par.range[0] = this.minv
    this.parameter_value_change()
  }

  max_change () {
    var val = this.par.val / this.par.factor;
    var max = this.maxv;
    if (isNaN(max)) {
      alert("Invalid input");
      this.change_parameter(this.par);
      return;
    }
    if (this.par.bounds) {
      if (!(max >= this.par.bounds[0] && max <= this.par.bounds[1])) {
        alert(`Input out of range: [${this.par.bounds}]`);
        this.change_parameter(this.par);
        return;
      }
    }
    if (this.maxv < this.minv) {
      this.min.val(this.max.val());
    }
    if (val > this.maxv) {
      this.par.store(this.maxv)
      this.par.set_field();
      var pos = this.pos_from_par(this.maxv);
    }
    else {
      var pos = this.pos_from_par()
    }
    this.par.range[1] = this.maxv
    this.parameter_value_change()
  }
}


function initAll() {
  slider = new Slider("#slider-bar", "#slider-min", "#slider-max");
  spins = initSpins(NUMBER_SPINS);
  noise = new Parameter("#experimentalNoise", "acquisition", "ufloat", 1.0, [0,1]);
  acquP = new AcquisitionParameters();
  procP = new ProcessParameters();

  tdPlot = initPlot("#plot-time-domain-svg");
  fqPlot = initPlot("#plot-freq-domain-svg");
  vcPlot = initVectorPlot()
  initMouseOverEffects(tdPlot, vcPlot)

  update();
}


function update(parameter_type) {
  initSignalArray(acquP.numberPoints);
  initTimeArray(procP.zeroFilling.val);
  initProcessArray();
  calculateSignal();
  processTimeDomain();
  updateTimeDomainPlot(tdPlot);
  processFrequencyDomain();
  updateFreqDomainPlot(fqPlot);
  updateProcessingPlots(tdPlot, fqPlot);
}


class ProcessParameters {
  constructor () {
    this.zeroFilling = new Parameter("#zeroFilling", "processing", "uint", 1.0, [2,2048] ,[0, MAX_ZEROFILL])
    this.phcor0 = new Parameter("#phaseCorrection0", "processing", "float", PI/180, [-180,180])
    this.phcor1 = new Parameter("#phaseCorrection1", "processing", "float", PI/180, [-360,360])

    this.windows = {
      none: {},
      exponential: {
        a: new Parameter("#exp_a", "processing", "ufloat", 1.0, [0,20])
      },
      gaussian: {
        a: new Parameter("#gauss_a", "processing", "ufloat", 1.0, [0,20]),
      },
      trigonometric: {
        a: new Parameter("#trig_a", "processing", "ufloat", 1.0, [0,2]),
        b: new Parameter("#trig_b", "processing", "ufloat", 1.0, [0,1]),
        c: new Parameter("#trig_c", "processing", "uint", 1.0, [0,4]),
      }
    }

    this.windowName = $("#windowFunction").val()
    $("#windowFunction").on("change", e => {
      $(`.window.${this.windowName}`).hide(100);
      this.windowName = e.target.value;
      $(`.window.${this.windowName}`).show(100);
    })
  }

  get window () {
    switch (this.windowName) {
      case "none":
        return false;
      break;
      case "exponential":
        return this.exp_func();
      break;
      case "gaussian":
        return this.gauss_func();
      break;
      case "trigonometric":
        return this.trig_func();
      break;
    }
  }

  exp_func () {
    let a = this.windows.exponential.a.val;
    return function (t) {
      return Math.exp(-1 * PI * a * t)
    }
  }

  gauss_func () {
    let a = this.windows.gaussian.a.val;
    return function (t) {
      return Math.exp(-1*(PI * a * t)**2)
    }
  }

  trig_func () {
    let a = this.windows.trigonometric.a.val;
    let b = this.windows.trigonometric.b.val;
    let c = this.windows.trigonometric.c.val;
    let taq = acquP.dwellTime * sig.length;
    return function (t) {
      return Math.sin(PI*a + PI*(b-a) * (t/taq))**c
    }
  }
}


class AcquisitionParameters {
  constructor () {
    this.variables = ["numberPoints","acquisitionTime","spectralWidth","dwellTime"];

    this.pars = {
      numberPoints: new Parameter("#numberPoints", "acquisition", "uint", 1.0, [2,2048], [2, MAX_POINTS]),
      acquisitionTime: new Parameter("#acquisitionTime", "acquisition", "ufloat", 1.0, [0,1]),
      dwellTime: new Parameter("#dwellTime", "acquisition", "ufloat", 1.0, [0.0001, 0.1]),
      spectralWidth: new Parameter("#spectralWidth", "acquisition", "ufloat", 1.0, [10,500]),
      carrierFrequency: new Parameter("#carrierFrequency", "acquisition", "float", 1.0, [-50,50]),
    }

    this.variables.forEach( v => {
      let par = this.pars[v];
      par.on_store = () => {this.setVariable(v)};
    })
  }

  get numberPoints () {
    return this.pars.numberPoints.val;
  }
  get acquisitionTime () {
    return this.numberPoints * this.dwellTime;
  }
  get dwellTime () {
    return this.pars.dwellTime.val;
  }
  get spectralWidth () {
    return 1.0 / this.dwellTime;
  }
  get carrierFrequency () {
    return this.pars.carrierFrequency.val;
  }

  quadrature () {
    return $("#quadratureDetection").is($(":checked"));
  }

  setVariable (variable) {
    switch (variable) {
      case "numberPoints":
      this.pars.acquisitionTime.store(this.acquisitionTime, false);
      this.pars.acquisitionTime.set_field();
      break;
      case "acquisitionTime":
      this.pars.numberPoints.store(this.pars.acquisitionTime.val / this.dwellTime, false);
      this.pars.numberPoints.set_field();
      break;
      case "dwellTime":
      this.pars.spectralWidth.store(this.spectralWidth, false);
      this.pars.spectralWidth.set_field();
      this.pars.acquisitionTime.store(this.acquisitionTime, false);
      this.pars.acquisitionTime.set_field();
      break;
      case "spectralWidth":
      this.pars.dwellTime.store(1.0 / this.pars.spectralWidth.val, false);
      this.pars.dwellTime.set_field();
      this.pars.acquisitionTime.store(this.acquisitionTime, false);
      this.pars.acquisitionTime.set_field();
      break;
    }
  }
}


function calculateSignal () {
  if (noise.val) {
    let scale = noise.val;
    for (let i=0; i<sig.length; i++) {
      sig.real[i] = scale*(Math.random()-0.5);
      sig.imag[i] = scale*(Math.random()-0.5);
    }
  }
  else {
    sig.real.fill(0);
    sig.imag.fill(0);
  }

  spins.forEach( spin => {
    let amp = spin.amplitude.val;
    let offs = spin.offset.val - acquP.carrierFrequency*2*PI;
    let ph = spin.phase.val;
    let t2 = spin.t2.val;
    let dt = acquP.dwellTime;
    for (let i=0; i<sig.length; i++) {
      let pf = amp * Math.exp(-(i*dt)/t2)
      let arg = (offs * i * dt) + ph
      sig.real[i] += pf * Math.cos(arg)
      sig.imag[i] += pf * Math.sin(arg)
    }
  })
}


function processTimeDomain () {

  proc.real.set(time.real);
  if (acquP.quadrature()) {
    proc.imag.set(time.imag);
  }
  else {
    proc.imag.fill(0.0);
  }

  let wfunc = procP.window
  if (wfunc) {
    let dt = acquP.dwellTime
    let func = procP.window;
    for (let i=0; i<sig.length; i++) {
      let a = wfunc(i * dt);
      proc.real[i] *= a;
      proc.imag[i] *= a;
    }
  }
}


function processFrequencyDomain () {

  proc.FFT();
  proc.fftShift();

  let ph0 = procP.phcor0.val;
  let ph1 = procP.phcor1.val;
  let shift = 0.5 * (proc.length % 2);

  for (let i=0; i<proc.length; i++) {
    let f = (0.5 - (i + shift)/proc.length);
    let theta = -(ph0 + f * ph1)
    let re = proc.real[i];
    let im = proc.imag[i];
    let cos = Math.cos(theta);
    let sin = Math.sin(theta);
    proc.real[i] = cos * re + sin * im;
    proc.imag[i] = cos * im - sin * re;
  }
}


function backCalculateSignal () {
  proc.fftShift();
  proc.InvFFT();
}


function initPlot (canvas_id) {

  var canvas = d3.select(canvas_id);
  var svg = canvas.append("g").attr("class", "frame")
  var margin = {top: 20, right: 20, bottom: 40, left: 30};
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
    .attr("class", "line_processing")

  svg.append("path")
    .attr("class", "line_imag")

  svg.append("path")
    .attr("class", "line_real")

  if (canvas_id=="#plot-time-domain-svg") {
    var xLabel = "Time /s";
  }
  else {
    var xLabel = "Frequency /Hz";
  }

  svg.append("text")
    .attr("y", height + 35)
    .attr("x", width / 2)
    .style("text-anchor", "middle")
    .style("font-size", "10pt")
    .text(xLabel); 

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
    proc: svg.select(".line_processing"),
    xaxis: svg.select(".x_axis"),
    yaxis: svg.select(".y_axis"),
    range: null
  }

  return plotData
}


function updateProcessingPlots (tdPlot, fqPlot) {

  let wfunc = procP.window;
  if (wfunc) {
    var wfdata = new Float64Array(PROC_PLOT_POINTS);
    let taq = acquP.dwellTime * sig.length;
    for (let i=0; i<PROC_PLOT_POINTS; i++) {
      let t = (i/PROC_PLOT_POINTS) * taq
      wfdata[i] = wfunc(t);
    }
  } else {
    var wfdata = [0];
  }

  let fac = (acquP.dwellTime * sig.length) / wfdata.length;
  var range = tdPlot.range;
  let wline = d3.line()
    .x( (d, i) => tdPlot.xScale(fac * i) )
    .y(  d     => tdPlot.yScale(d * range) )

  tdPlot.proc.datum(wfdata)
    .transition()
    .duration(ANI_DUR)
    .attr("d", wline)

  let ph0 = procP.phcor0.val;
  let ph1 = procP.phcor1.val;
  let pi2 = 2*PI
  if (ph0 || ph1) {
    var phdata = new Float64Array(PROC_PLOT_POINTS);
    let sw = acquP.spectralWidth;
    for (let i=0; i<PROC_PLOT_POINTS; i++) {
      let f = 0.5 - i/PROC_PLOT_POINTS;
      phdata[i] = (((((ph0 + f * ph1) / PI) - 1) % 2) + 2) % 2 - 1;
    }
  }
  else {
    var phdata = [0];
  }
  
  let sw = acquP.spectralWidth;
  var range = fqPlot.range;
  let pline = d3.line()
    .x( (d, i) => fqPlot.xScale((0.5 - i/phdata.length)*sw + acquP.carrierFrequency))
    .y(  d     => fqPlot.yScale(d * range) )

  fqPlot.proc.datum(phdata)
    .transition()
    .duration(ANI_DUR)
    .attr("d", pline)
}


function updateTimeDomainPlot (plot) {

  let maximini = [
    Math.min.apply(null, proc.real),
    Math.max.apply(null, proc.real),
    Math.min.apply(null, proc.imag),
    Math.max.apply(null, proc.imag),
  ];
  let range = Math.max(...maximini.map( v => Math.abs(v)));
  plot.range = range;
  plot.yScale.domain([-range, range])

  plot.yaxis.transition()
    .duration(ANI_DUR)
    .call(d3.axisLeft(plot.yScale))

  plot.xScale.domain([0, acquP.dwellTime * proc.length])
  plot.xaxis.transition()
    .duration(ANI_DUR)
    .call(d3.axisBottom(plot.xScale))

  const dt = acquP.dwellTime
  let line = d3.line()
    .x( (d, i) => plot.xScale(i*dt) )
    .y(  d     => plot.yScale(d) );

  plot.real.datum(proc.real)
    .transition()
    .duration(ANI_DUR)
    .attr("d", line)

  plot.imag.datum(proc.imag)
    .transition()
    .duration(ANI_DUR)
    .attr("d", line)
}


function updateFreqDomainPlot (plot) {

  let maximini = [
    Math.min.apply(null, proc.real),
    Math.max.apply(null, proc.real),
    Math.min.apply(null, proc.imag),
    Math.max.apply(null, proc.imag),
  ];
  let range = Math.max(...maximini.map( v => Math.abs(v)))
  plot.range = range;
  plot.yScale.domain([-range, range])
  plot.yaxis.transition()
    .duration(ANI_DUR)
    .call(d3.axisLeft(plot.yScale))

  plot.xScale.domain([0.5*acquP.spectralWidth + acquP.carrierFrequency,
                     -0.5*acquP.spectralWidth + acquP.carrierFrequency])
  plot.xaxis.transition()
    .duration(ANI_DUR)
    .call(d3.axisBottom(plot.xScale))

  let shift = 0.5 * (proc.length % 2);
  let line = d3.line()
    .x( (d, i) => plot.xScale(
      (0.5 - (i + shift)/proc.length) * acquP.spectralWidth + acquP.carrierFrequency))
    .y(  d     => plot.yScale(d) );

  plot.real.datum(proc.real)
    .transition()
    .duration(ANI_DUR)
    .attr("d", line)

  plot.imag.datum(proc.imag)
    .transition()
    .duration(ANI_DUR)
    .attr("d", line)
}


function initVectorPlot () {

  var canvas = d3.select("#plot-spin-vector-svg");
  var svg = canvas.append("g").attr("class", "frame")
  var margin = {top: 40, right: 40, bottom: 20, left: 20};
  var width = canvas.attr("width") - margin.left - margin.right;
  var height = canvas.attr("height") - margin.top - margin.bottom;
  var radius = 0.5*height;
  
  svg.attr("transform", `translate(${margin.left+0.5*width}, ${margin.top+0.5*height})`);

  svg.append("circle")
    .attr("class", "circle")
    .attr("r", radius)
    .attr("fill", "none");

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
      .attr("class", `td-mouse-over vector spin spin${i}`)
      .attr("visibility", "hidden")
    spinG.append("line")
      .attr("marker-end", appendArrowheadMarker(svg, "spin"))

    let pointerG = spinG.append("g")
      .attr("class", `pointer`)
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

  let legh = 20;
  let legw = 70;

  svg.append("g")
    .attr("class", "td-mouse-over vector sum")
    .attr("visibility", "hidden")
    .append("line")
    .attr("marker-end", appendArrowheadMarker(svg, "sum"));

  let realG = svg.append("g")
    .attr("transform", `translate(80, ${-0.5*legh-100})`)
    .attr("class", "legend real")

  realG.append("text")
    .text("real:")
    .attr("fill", "black")
    .attr("font-size", "14")
    .attr("alignment-baseline", "middle")

  realG.append("line")
    .attr("class", "line_real")
    .attr("x1", 35)
    .attr("x2", legw)

  realG.append("rect")
    .attr("fill", "none")
    .attr("opacity", 0.5)
    .attr("pointer-events", "all")
    .attr("width", legw)
    .attr("height", legh)
    .attr("y", -legh*0.5)
    .on("click", () => {
      let sele = d3.selectAll(".line_real")
      let opac = parseFloat(sele.style("opacity"))
      if (opac < 0.5){
         sele.transition()
          .duration(ANI_DUR)
          .style("opacity", 1);
      }
      else {
        sele.transition()
          .duration(ANI_DUR)
          .style("opacity", 0.15);
      }
    })


  let imagG = svg.append("g")
    .attr("transform", `translate(80, ${0.5*legh-100})`)
    .attr("class", "legend imag")

  imagG.append("text")
    .text("imag:")
    .attr("fill", "black")
    .attr("font-size", "14")
    .attr("alignment-baseline", "middle")

  imagG.append("line")
    .attr("class", "line_imag")
    .attr("x1", 35)
    .attr("x2", legw)

  imagG.append("rect")
    .attr("fill", "none")
    .attr("opacity", 0.5)
    .attr("pointer-events", "all")
    .attr("width", legw)
    .attr("height", legh)
    .attr("y", -legh*0.5)
    .on("click", () => {
      let sele = d3.selectAll(".line_imag")
      let opac = parseFloat(sele.style("opacity"))
      if (opac < 0.5){
         sele.transition()
          .duration(ANI_DUR)
          .style("opacity", 1);
      }
      else {
        sele.transition()
          .duration(ANI_DUR)
          .style("opacity", 0.15);
      }
    })



  return {
    svg: svg,
    radius: radius
  }
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


function initMouseOverEffects (tdPlot, vcPlot) {
  var mouseG = tdPlot.svg.append("g")
    .attr("class", "td-mouse-over-effects");

  mouseG.append("path")
    .attr("class", "td-mouse-over mouse-line")
    .attr("visibility", "hidden")

  mouseG.append('svg:rect') // append a rect to catch mouse movements on canvas
    .attr('width', tdPlot.width) // can't catch mouse events on a g element
    .attr('height', tdPlot.height)
    .attr('fill', 'none')
    .attr('pointer-events', 'all')
    .on('mouseout', () => {
      d3.selectAll(".td-mouse-over")
        .attr("visibility", "hidden");
    })
    .on('mouseover', () => {
      d3.selectAll(".td-mouse-over")
        .attr("visibility", "visible");
      spins.forEach( (spin, i) => {
        if (!spin.amplitude.val) {
          d3.selectAll(`#plot-spin-vector-svg .spin${i}`)
            .attr("visibility", "hidden")
        }
      })
    })
    .on('mousemove', () => {
      let mouse = d3.mouse(d3.event.currentTarget);
      let t_abs = tdPlot.xScale.invert(mouse[0])
      let p = Math.abs(Math.round( t_abs / acquP.dwellTime))
      let t = p * acquP.dwellTime;
      let loc = tdPlot.xScale(t);
      if (p < proc.length) {
        d3.select(".mouse-line")
          .attr("d", () => `M${loc},${tdPlot.height} ${loc},0`);
        d3.select("#plot-spin-vector-svg .vector.sum")
          .select("line")
          .attr("x2", (time.real[p] / tdPlot.range) * vcPlot.radius)
          .attr("y2", -(time.imag[p] / tdPlot.range) * vcPlot.radius)
        spins.forEach( (spin, i) => {
          let mag = spin.evolution(t_abs);
          let theta = Math.atan2(mag.real, mag.imag) - PI/2
          d3.select(`#plot-spin-vector-svg .vector.spin${i}`)
            .select("line")
            .attr("x2", (mag.real / tdPlot.range) * vcPlot.radius)
            .attr("y2", -(mag.imag / tdPlot.range) * vcPlot.radius);
          d3.select(`#plot-spin-vector-svg .vector.spin${i} .pointer`)
            .attr("transform", `rotate(${(180.0/PI)*theta})`)
        })
      }
    })
}






































