
// Global variables
const ArrayType = Float64Array;
const MAX_POINTS = (2**12) // 16k data
const MAX_ZEROFILL = (2**12) // 16k zerofill
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

  deregister () {
    this.elem.off();
  }

  store(value=this.elem.val()) {
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
    return true
  }

  set_field(value=this.val) {
    this.elem.val(value / this.factor);
  }

  set_value(value) {
    this.val = value;
    this.set_field();
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
    let id = this.par.id.slice(1)
    if (acquP.variables.includes(id)) {
      acquP.setVariable(id);
    }
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
  initTimeArray(procP.zeroFilling.val)
  initProcessArray()
  calculateSignal();
  processSignal();
  updateFreqDomainPlot(fqPlot);
  proc.fftShift();
  proc.InvFFT();
  updateTimeDomainPlot(tdPlot);
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
      gaussian: {},
      trigonometric: {}
    }

    this.windowName = $("#windowFunction").val()
    $("#windowFunction").on("change", e => {
      $(`.window.${this.windowName}`).hide(1000);
      this.windowName = e.target.value;
      $(`.window.${this.windowName}`).show(1000);
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
    }
  }

  exp_func () {
    let a = this.windows.exponential.a.val
    return function (t) {
      return Math.exp(-PI * a * t)
    }
  }

  gauss_func () {
    let a = this.windows.gaussian.a.val
    return function (t) {
      return Math.exp( -a*PI*t - b*t*t)
    }
  }



  // EM[i] = exp( -PI*i*lb/sw )

  // SP[i] = sin( PI*off + PI*(end-off)*i/(tSize-1) )^pow

//           GMB[i] = exp( -PI*lb*t + (PI*lb/(2.0*gb*aq))*t*t )
//           t  = i/sw
//           aq = tSize/sw
//           a  = PI*lb
//           b  = -a/(2.0*gb*aq)



// Given the GMB adjustable parameters lb and gb:

//           GMB[i] = exp( -a*t - b*t*t )
// where
//           t  = i/sw
//           aq = tSize/sw
//           a  = PI*lb
//           b  = -a/(2.0*gb*aq)

}

class AcquisitionParameters {
  constructor () {
    this.lastSet = "numberPoints";
    this.lastLastSet = "acquisitionTime";
    this.variables = ["numberPoints","acquisitionTime","spectralWidth","dwellTime"];

    this.pars = {
      numberPoints: new Parameter("#numberPoints", "acquisition", "uint", 1.0, [2,2048], [0, MAX_POINTS]),
      acquisitionTime: new Parameter("#acquisitionTime", "acquisition", "ufloat", 1.0, [0,1]),
      dwellTime: new Parameter("#dwellTime", "acquisition", "ufloat", 1.0, [0.0001, 0.1]),
      spectralWidth: new Parameter("#spectralWidth", "acquisition", "ufloat", 1.0, [10,500]),
      carrierFrequency: new Parameter("#carrierFrequency", "acquisition", "float", 1.0, [-50,50]),
    }
    this.numberPoints = this.pars.numberPoints.val
    this.acquisitionTime = this.pars.acquisitionTime.val

    this.variables.forEach( v => {
      let par = this.pars[v]
      par.elem.off();
      par.elem.on("keypress", e => {
        if (e.which==13) {
          if (par.store()) {
            this.setVariable(v);
            update(par.parameter_type);
          }
        };
      })
      this.pars[v].elem.on("focusout", e => {
        this.setVariable(v);
      })
    })
  }
  get dwellTime () {
    return this.acquisitionTime / this.numberPoints;
  }
  get spectralWidth () {
    return 1.0 / this.dwellTime;
  }
  get carrierFrequency () {
    return this.pars.carrierFrequency.val;
  }
  set carrierFrequency (value) {
    this.pars.carrierFrequency.val = value;
  }

  quadrature () {
    return $("#quadratureDetection").is($(":checked"));
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
    this.pars.numberPoints.set_field(this.numberPoints);
    this.pars.acquisitionTime.set_field(this.acquisitionTime);
    this.pars.dwellTime.set_field(this.dwellTime);
    this.pars.spectralWidth.set_field(this.spectralWidth)
  }
  setVariable (variable) {
    switch (variable) {
      case "numberPoints":
      this.setNumberPoints(this.pars.numberPoints.val);
      break;
      case "acquisitionTime":
      this.setAcquisitionTime(this.pars.acquisitionTime.val);
      break;
      case "dwellTime":
      this.setDwellTime(this.pars.dwellTime.val);
      break;
      case "spectralWidth":
      this.setSpectralWidth(this.pars.spectralWidth.val);
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
      this.numberPoints = math.round(value / this.dwellTime);
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

function processSignal () {
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
    .attr("class", "line_processing")

  svg.append("path")
    .attr("class", "line_imag")

  svg.append("path")
    .attr("class", "line_real")

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
    let taq = acquP.dwellTime * proc.length;
    for (let i=0; i<PROC_PLOT_POINTS; i++) {
      let t = (i/PROC_PLOT_POINTS) * taq
      wfdata[i] = wfunc(t);
    }
  } else {
    var wfdata = [0];
  }

  let fac = (acquP.dwellTime * proc.length) / wfdata.length;
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

  svg.append("g")
    .attr("class", "td-mouse-over vector sum")
    .attr("visibility", "hidden")
    .append("line")
    .attr("marker-end", appendArrowheadMarker(svg, "sum"));

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
      let p = Math.abs(Math.round( tdPlot.xScale.invert(mouse[0]) / acquP.dwellTime))
      let t = p * acquP.dwellTime;
      let loc = tdPlot.xScale(t);
      if (p < proc.length) {
        d3.select(".mouse-line")
          .attr("d", () => `M${loc},${tdPlot.height} ${loc},0`);
        d3.select("#plot-spin-vector-svg .vector.sum")
          .select("line")
          .attr("x2", (proc.real[p] / tdPlot.range) * vcPlot.radius)
          .attr("y2", -(proc.imag[p] / tdPlot.range) * vcPlot.radius)
        spins.forEach( (spin, i) => {
          let mag = spin.evolution(t);
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










// function updateTimeDomainPlot (plot) {

//   let range = math.max([math.abs(math.re(timeDomain)), math.abs(math.im(timeDomain))])
//   plot.range = range;
//   plot.yScale.domain([-range, range])
//   plot.yaxis.transition()
//     .duration(animationDuration)
//     .call(d3.axisLeft(plot.yScale))

//   plot.xScale.domain([0, acquisitionTime])
//   plot.xaxis.transition()
//     .duration(animationDuration)
//     .call(d3.axisBottom(plot.xScale))

//   let line = d3.line()
//     .x( (d, i) => plot.xScale(timeAxis[i]) )
//     .y(  d     => plot.yScale(d) );

//   plot.real.datum(math.re(timeDomain))
//     .transition()
//     .duration(animationDuration)
//     .attr("d", line)

//   plot.imag.datum(math.im(timeDomain))
//     .transition()
//     .duration(animationDuration)
//     .attr("d", line)

// }


// function updateFreqDomainPlot (plot) {

//   plot.yScale.domain(d3.extent(math.re(freqDomain).concat(math.im(freqDomain))))
//   plot.yaxis.transition()
//     .duration(animationDuration)
//     .call(d3.axisLeft(plot.yScale))

//   plot.xScale.domain([0.5*spectralWidth + carrierFrequency, -0.5*spectralWidth + carrierFrequency])
//   plot.xaxis.transition()
//     .duration(animationDuration)
//     .call(d3.axisBottom(plot.xScale))

//   let line = d3.line()
//     .x( (d, i) => plot.xScale(freqAxis[i]) )
//     .y(  d     => plot.yScale(d) );

//   plot.real.datum(math.re(freqDomain))
//     .transition()
//     .duration(animationDuration)
//     .attr("d", line)

//   plot.imag.datum(math.im(freqDomain))
//     .transition()
//     .duration(animationDuration)
//     .attr("d", line)

// }

// function update () {
//   t0 = performance.now()
//   aPars.storeCurrent();
//   populateVariables();
//   initAxes();
//   calculateSignal();
//   processSignal();
//   updateTimeDomainPlot(tdPlot);
//   updateFreqDomainPlot(fqPlot);
//   t1 = performance.now()
//   console.log(t1-t0)
// }


// (function($) {
//   $.fn.inputFilter = function(inputFilter) {
//     return this.on("input keydown keyup mousedown mouseup select contextmenu drop", function() {
//       if (inputFilter(this.value)) {
//         this.oldValue = this.value;
//         this.oldSelectionStart = this.selectionStart;
//         this.oldSelectionEnd = this.selectionEnd;
//       } else if (this.hasOwnProperty("oldValue")) {
//         this.value = this.oldValue;
//         this.setSelectionRange(this.oldSelectionStart, this.oldSelectionEnd);
//       } else {
//         this.value = "";
//       }
//     });
//   };
// }(jQuery));


// function initInputFilters () {
//   $("#numberPoints").inputFilter( value => /^\d*$/.test(value) );
//   $("#acquisitionTime").inputFilter( value => /^\d*[.]?\d*$/.test(value) );
//   $("#spectralWidth").inputFilter( value => /^\d*[.]?\d*$/.test(value) );
//   $("#dwellTime").inputFilter( value => /^\d*[.]?\d*$/.test(value) );
//   $("#carrierFrequency").inputFilter( value => /^-?\d*[.]?\d*$/.test(value) );

//   $(".spins .amplitude").each( (i, elem) => {
//     $(elem).inputFilter( value => /^\d*[.]?\d*$/.test(value) )
//   })
//   $(".spins .offset").each( (i, elem) => {
//     $(elem).inputFilter( value => /^-?\d*[.]?\d*$/.test(value) )
//   })
//   $(".spins .t2").each( (i, elem) => {
//     $(elem).inputFilter( value => /^\d*[.]?\d*$/.test(value) )
//   })
//   $(".spins .phase").each( (i, elem) => {
//     $(elem).inputFilter( value => /^-?\d*[.]?\d*$/.test(value) )
//   })
// }


// $(document).ready( () => {
//   $(document).on('keypress',function(e) {
//     if(e.which == 13) {
//       update();
//     }
//   });

//   initInputFilters();
//   update();

 

// });


































