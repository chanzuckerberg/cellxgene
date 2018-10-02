// jshint esversion: 6
const createCamera = require("orbit-camera");
const createScroll = require("scroll-speed");
const mp = require("mouse-position");
const mb = require("mouse-pressed");
const key = require("key-pressed");

const panSpeed = 0.4;
const scaleSpeed = 0.5;
const scaleMax = 3;
// const scaleMin = 1.15
const scaleMin = 1.03;

function attachCamera(canvas, opts) {
  opts = opts || {};
  opts.pan = opts.pan !== false;
  opts.scale = opts.scale !== false;
  opts.rotate = opts.rotate !== false;

  const scroll = createScroll(canvas, opts.scale);
  const mbut = mb(canvas, opts.rotate);
  const mpos = mp(canvas);
  const camera = createCamera([0, 0, 1], [0, 0, -1], [0, 1, 0]);

  camera.tick = tick;

  return camera;

  function tick() {
    const ctrl = key("<control>") || key("<alt>");
    const alt = key("<shift>");
    const { height, width } = canvas;

    if (opts.rotate && mbut.left && ctrl && !alt) {
      camera.rotate(
        [mpos.x / width - 0.5, mpos.y / height - 0.5],
        [mpos.prevX / width - 0.5, mpos.prevY / height - 0.5]
      );
    }

    if ((opts.pan && mbut.right) || (mbut.left && !ctrl && !alt)) {
      camera.pan([
        ((panSpeed * (mpos[0] - mpos.prev[0])) / width) * camera.distance,
        ((panSpeed * (mpos[1] - mpos.prev[1])) / height) * camera.distance
      ]);
    }

    if (opts.scale && scroll[1]) {
      camera.distance *= Math.exp((scroll[1] * scaleSpeed) / height);
    }

    if (opts.scale && (mbut.middle || (mbut.left && !ctrl && alt))) {
      const d = mpos.y - mpos.prevY;
      if (!d) return;

      camera.distance *= Math.exp(d / height);
    }

    if (camera.distance > scaleMax) camera.distance = scaleMax;
    if (camera.distance < scaleMin) camera.distance = scaleMin;

    scroll.flush();
    mpos.flush();
  }
}

export default attachCamera;
