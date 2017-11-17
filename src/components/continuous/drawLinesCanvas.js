import {
  project
} from "./util";

import renderQueue from "./renderQueue";

/*****************************************
******************************************
              draw loop
******************************************
******************************************/

const drawLinesCanvas = (ctx, dimensions, xscale) => {
  return (d) => {

    // ctx.strokeStyle = d["Sample.name.color"];
    ctx.strokeStyle = "rgba(0,0,0,1)";
    ctx.beginPath();
    var coords = project(d, dimensions, xscale);
    coords.forEach((p,i) => {
      // this tricky bit avoids rendering null values as 0
      if (p === null) {
        // this bit renders horizontal lines on the previous/next
        // dimensions, so that sandwiched null values are visible
        if (i > 0) {
          var prev = coords[i-1];
          if (prev !== null) {
            ctx.moveTo(prev[0],prev[1]);
            ctx.lineTo(prev[0]+6,prev[1]);
          }
        }
        if (i < coords.length-1) {
          var next = coords[i+1];
          if (next !== null) {
            ctx.moveTo(next[0]-6,next[1]);
          }
        }
        return;
      }

      if (i == 0) {
        ctx.moveTo(p[0],p[1]);
        return;
      }

      ctx.lineTo(p[0],p[1]);
    });
    ctx.stroke();
  }
}

const drawCellLinesUsingRenderQueue = (
  metadata,
  dimensions,
  xscale,
  ctx,
  width,
  height
) => {

  const _renderLinesWithQueue = renderQueue(
    drawLinesCanvas(ctx, dimensions, xscale)
  ).rate(50);
  
  _renderLinesWithQueue(metadata);

  return _renderLinesWithQueue;

}

export default drawCellLinesUsingRenderQueue;
