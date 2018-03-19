var createCamera = require('orbit-camera')
var createScroll = require('scroll-speed')
var mp = require('mouse-position')
var mb = require('mouse-pressed')
var key = require('key-pressed')

const panSpeed = 0.4
const scaleSpeed = 0.5
const scaleMax = 3
const scaleMin = 1.15

function attachCamera(canvas, opts) {
  opts = opts || {}
  opts.pan = opts.pan !== false
  opts.scale = opts.scale !== false
  opts.rotate = opts.rotate !== false

  var scroll = createScroll(canvas, opts.scale)
  var mbut = mb(canvas, opts.rotate)
  var mpos = mp(canvas)
  var camera = createCamera(
      [0, 0, 1]
    , [0, 0, -1]
    , [0, 1, 0]
  )

  camera.tick = tick

  return camera

  function tick() {
    var ctrl = key('<control>') || key('<alt>')
    var alt = key('<shift>')
    var height = canvas.height
    var width = canvas.width

    if (opts.rotate && mbut.left && ctrl && !alt) {
      camera.rotate(
          [ mpos.x / width - 0.5, mpos.y / height - 0.5 ]
        , [ mpos.prevX / width - 0.5, mpos.prevY / height - 0.5 ]
      )
    }

    if (opts.pan && mbut.right || (mbut.left && !ctrl && !alt)) {
      camera.pan([
          (panSpeed * (mpos[0] - mpos.prev[0]) / width) * Math.pow(camera.distance, 1)
        , (panSpeed * (mpos[1] - mpos.prev[1]) / height) * Math.pow(camera.distance, 1)
      ])
    }

    if (opts.scale && scroll[1]) {
      camera.distance *= Math.exp(scroll[1] * scaleSpeed / height)
    }

    if (opts.scale && (mbut.middle || (mbut.left && !ctrl && alt))) {
      var d = (mpos.y - mpos.prevY)
      if (!d) return;

      camera.distance *= Math.exp(d / height)
    }

    if (camera.distance > scaleMax) camera.distance = scaleMax
    if (camera.distance < scaleMin) camera.distance = scaleMin

    scroll.flush()
    mpos.flush()
  }
}

export default attachCamera;
