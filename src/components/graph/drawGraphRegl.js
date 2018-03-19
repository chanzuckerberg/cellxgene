import mat4 from 'gl-mat4';
import fit from 'canvas-fit';
import _camera from '../../util/camera.js'
import _regl from 'regl'
import _drawPoints from './drawPointsRegl'

// setup canvas and camera
const canvas = document.body.appendChild(document.createElement('canvas'))
console.log('canvas',canvas)
// css(canvas, {zIndex: -1000})

const camera = _camera(canvas, {scale: true, rotate: false});
const regl = _regl(canvas)

window.addEventListener('resize', fit(canvas), false)

// add explanation
const div = document.createElement('h1')
// div.innerHTML = 'press spacebar to change'
// css(div, {'font-family': 'monospace', 'text-align': 'center'})
// document.body.appendChild(div)

// import draw functions
const drawPoints = _drawPoints(regl)

// import data generators
const generatePoints = function (count) {
  return Array(count).fill().map(function () {
    return [
      Math.random() - 0.5,
      Math.random() - 0.5
    ]
  })
}

const generateColors = function (count) {
  return Array(count).fill().map(function () {
    return [
      Math.random(),
      Math.random(),
      Math.random()
    ]
  })
}

// set constants
const count = 3e4

// preallocate buffers
const pointBuffer = regl.buffer(count)(generatePoints(count))
const colorBuffer = regl.buffer(count)(generateColors(count))

// update on spacebar
document.body.onkeyup = function(e){
  if (e.keyCode == 32){
    pointBuffer(generatePoints(count))
    colorBuffer(generateColors(count))
  }
}

regl.frame(({viewportWidth, viewportHeight, time}) => {

  regl.clear({
    depth: 1,
    color: [1, 1, 1, 1]
  })

  drawPoints({
    distance: camera.distance,
    color: colorBuffer,
    position: pointBuffer,
    count: count,
    view: camera.view(),
    scale: viewportHeight/viewportWidth
  })

  camera.tick()
})
