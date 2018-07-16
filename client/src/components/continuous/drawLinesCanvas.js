// jshint esversion: 6
import _ from "lodash";
import { project } from "./util";

import renderQueue from "../../util/renderQueue";

/*****************************************
******************************************
              draw loop
******************************************
******************************************/

const drawLinesCanvas = (
  ctx,
  dimensions,
  xscale,
  colorAccessor,
  colorScale
) => {
  return d => {
    ctx.globalAlpha = 0.1;

    if (d["__selected__"]) {
      ctx.strokeStyle = d["__color__"];
    } else {
      return;
    }

    ctx.beginPath();
    var coords = project(d, dimensions, xscale);
    coords.forEach((p, i) => {
      // this tricky bit avoids rendering null values as 0
      if (p === null) {
        // this bit renders horizontal lines on the previous/next
        // dimensions, so that sandwiched null values are visible
        if (i > 0) {
          var prev = coords[i - 1];
          if (prev !== null) {
            ctx.moveTo(prev[0], prev[1]);
            ctx.lineTo(prev[0] + 6, prev[1]);
          }
        }
        if (i < coords.length - 1) {
          var next = coords[i + 1];
          if (next !== null) {
            ctx.moveTo(next[0] - 6, next[1]);
          }
        }
        return;
      }

      if (i == 0) {
        ctx.moveTo(p[0], p[1]);
        return;
      }

      ctx.lineTo(p[0], p[1]);
    });
    ctx.stroke();
  };
};

const drawCellLinesUsingRenderQueue = (
  metadata,
  dimensions,
  xscale,
  ctx,
  colorAccessor,
  colorScale
) => {
  const _renderLinesWithQueue = renderQueue(
    drawLinesCanvas(ctx, dimensions, xscale, colorAccessor, colorScale)
  ).rate(50);
  _renderLinesWithQueue(metadata);
  return _renderLinesWithQueue;
};

const drawCellLinesSync = (
  metadata,
  dimensions,
  xscale,
  ctx,
  colorAccessor,
  colorScale
) => {
  const _draw = drawLinesCanvas(
    ctx,
    dimensions,
    xscale,
    colorAccessor,
    colorScale
  );
  _.each(metadata, _draw);
};

export default drawCellLinesUsingRenderQueue;
// export default drawCellLinesUsingRenderQueue;
