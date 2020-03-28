const math = require('mathjs');
const fft = require('jsfft');
const d3 = require('d3');


function fftShift(arr) {
  let shift = (arr.length + 1) / 2;
  shift -= arr.length * Math.floor(shift / arr.length);
  arr.push.apply(arr, arr.splice(0, shift));
  return arr;
}


class Spin {
  constructor(frequency=6*2*math.pi, r1=1.0, r2=1.0) {
    this.frequency = frequency;
    this.r1 = r1;
    this.r2 = r2;
  }
  exponent() {
    return math.complex(-this.r2, this.frequency);
  }
}

var spins = new Array(1).fill(new Spin())
var numberPoints;
var acquisitionTime;
var spectralWidth;
var timeAxis;
var freqAxis;
var timeDomain;
var freqDomain;
var fftdata;

function populateVariables () {
  numberPoints = parseFloat(document.getElementById("numberPoints").value);
  acquisitionTime = parseFloat(document.getElementById("acquisitionTime").value);
  spectralWidth = parseFloat(document.getElementById("spectralWidth").value);
}

function initAxes () {
  timeAxis = new Array(numberPoints);
  freqAxis = new Array(numberPoints);
  for (let i=0; i<numberPoints; i++) {
    timeAxis[i] = acquisitionTime * (i / numberPoints)
    freqAxis[i] = (0.5 - (i / numberPoints)) * spectralWidth
  }
}

function calculateSignal () {
  timeDomain = new Array(numberPoints).fill(math.complex(0.0));
  spins.forEach(
    spin => {
      timeDomain = math.add(timeDomain, math.exp(math.multiply(spin.exponent(), timeAxis)))
    }
  )
}

populateVariables();
initAxes();
calculateSignal();
drawPlot();




function drawPlot () {

  // 2. Use the margin convention practice 
  var margin = {top: 50, right: 50, bottom: 50, left: 50}
    , width = window.innerWidth - margin.left - margin.right // Use the window's width 
    , height = (window.innerHeight - margin.top - margin.bottom)*0.5; // Use the window's height

  // 5. X scale will use the index of our data
  var xScale = d3.scaleLinear()
      .domain([0, acquisitionTime])
      .range([0, width]);

  // 6. Y scale will use the randomly generate number 
  var yScale = d3.scaleLinear()
      .domain([-1, 1])
      .range([height, 0]);

  // 7. d3's line generator
  var line = d3.line()
      .x(function(d, i) { return xScale(timeAxis[i]); }) // set the x values for the line generator
      .y(function(d) { return yScale(d); }) // set the y values for the line generator 

  // 1. Add the SVG to the page and employ #2
  var svg = d3.select("body").append("svg")
      .attr("width", width + margin.left + margin.right)
      .attr("height", height + margin.top + margin.bottom)
    .append("g")
      .attr("transform", "translate(" + margin.left + "," + margin.top + ")");

  // 3. Call the x axis in a group tag
  svg.append("g")
      .attr("class", "x axis")
      .attr("transform", "translate(0," + height + ")")
      .call(d3.axisBottom(xScale)); // Create an axis component with d3.axisBottom

  // 4. Call the y axis in a group tag
  svg.append("g")
      .attr("class", "y axis")
      .call(d3.axisLeft(yScale)); // Create an axis component with d3.axisLeft

  // 9. Append the path, bind the data, and call the line generator 
  svg.append("path")
      .data([math.re(timeDomain)]) // 10. Binds data to the line 
      .attr("class", "line") // Assign a class for styling 
      .attr("d", line); // 11. Calls the line generator 


}
















// const c3 = require('c3');

// var chartTimeDomain = c3.generate({
//   bindto: '#chart-time-domain',
//   data: {
//     x: 'time',
//     columns: [
//       ['time'].concat(time),
//       ["Real"].concat(real), 
//       ["Imaginary"].concat(imag)
//     ]
//   },
//   axis: {
//     x: {
//       label: "Time /s",
//       tick: {
//         format: d3.format('.3n')
//       }
//     },
//     y: {
//       min: -1,
//       max: 1
//     }
//   },
//   point: {
//       show: false
//   },
// });

// var chartFreqDomain = c3.generate({
//   bindto: '#chart-freq-domain',
//   data: {
//     x: 'freq',
//     columns: [
//       ['freq'].concat(freq),
//       freal, 
//       fimag
//     ]
//   },
//   axis: {
//     x: {
//       label: "Time /s",
//       tick: {
//         format: d3.format('.3n')
//       }
//     },
//     // y: {
//     //   min: -1,
//     //   max: 1
//     // }
//   },
//   point: {
//       show: false
//   }
// });
