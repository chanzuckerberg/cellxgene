import * as d3 from "d3";
import Lasso from "./setupLasso";

/******************************************
*******************************************
          put svg & brush in DOM
*******************************************
******************************************/

// eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types --- FIXME: disabled temporarily on migrate to TS.
export default (
  // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any -- - FIXME: disabled temporarily on migrate to TS.
  selectionToolType: any,
  // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any -- - FIXME: disabled temporarily on migrate to TS.
  handleStartAction: any,
  // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any -- - FIXME: disabled temporarily on migrate to TS.
  handleDragAction: any,
  // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any -- - FIXME: disabled temporarily on migrate to TS.
  handleEndAction: any,
  // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any -- - FIXME: disabled temporarily on migrate to TS.
  handleCancelAction: any,
  // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any -- - FIXME: disabled temporarily on migrate to TS.
  viewport: any
) => {
  const svg = d3.select("#graph-wrapper").select("#lasso-layer");
  if (svg.empty()) return {};

  if (selectionToolType === "brush") {
    const brush = d3
      .brush()
      .extent([
        [0, 0],
        [viewport.width, viewport.height],
      ])
      .on("start", handleStartAction)
      .on("brush", handleDragAction)
      // FYI, brush doesn't generate cancel
      .on("end", handleEndAction);

    const brushContainer = svg
      .append("g")
      .attr("class", "graph_brush")
      .call(brush);

    return { svg, container: brushContainer, tool: brush };
  }

  if (selectionToolType === "lasso") {
    // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
    const lasso = (Lasso() as any)
      .on("end", handleEndAction)
      // FYI, Lasso doesn't generate drag
      .on("start", handleStartAction)
      .on("cancel", handleCancelAction);

    const lassoContainer = svg.call(lasso);

    return { svg, container: lassoContainer, tool: lasso };
  }

  throw new Error("unknown graph selection tool");
};
