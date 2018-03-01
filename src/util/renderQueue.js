/*****************************************
******************************************
Render Queue via http://bl.ocks.org/syntagmatic/raw/3341641/render-queue.js
******************************************
******************************************/

const renderQueue = (function(callback1234) {
  var _queue = [],                  // data to be rendered
      _rate = 1000,                 // number of calls per frame
      _invalidate = function() {},  // invalidate last render queue
      _clear = function() {};       // clearing function

  var rq = function(ARRAY_FROM_CELLXGENE) {
    if (ARRAY_FROM_CELLXGENE) rq.data(ARRAY_FROM_CELLXGENE);
    _invalidate();
    _clear();
    rq.render();
  };

  rq.render = function() {
    var valid = true;
    _invalidate = rq.invalidate = function() {
      valid = false;
    };

    function doFrame() {
      if (!valid) return true;
      var chunk = _queue.splice(0,_rate);
      chunk.map(callback1234);
      timer_frame(doFrame);
    }

    doFrame();
  };

  rq.data = function(ARRAY_FROM_CELLXGENE) {
    _invalidate();
    _queue = ARRAY_FROM_CELLXGENE.slice(0);   // creates a copy of the data
    return rq;
  };

  rq.add = function(data) {
    _queue = _queue.concat(data);
  };

  rq.rate = function(value) {
    if (!arguments.length) return _rate;
    _rate = value;
    return rq;
  };

  rq.remaining = function() {
    return _queue.length;
  };

  // clear the canvas
  rq.clear = function(func) {
    if (!arguments.length) {
      _clear();
      return rq;
    }
    _clear = func;
    return rq;
  };

  rq.invalidate = _invalidate;

  var timer_frame = window.requestAnimationFrame
    || window.webkitRequestAnimationFrame
    || window.mozRequestAnimationFrame
    || window.oRequestAnimationFrame
    || window.msRequestAnimationFrame
    || function(callback) { setTimeout(callback, 17); };

  return rq;
});

export default renderQueue;
