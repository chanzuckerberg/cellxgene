import * as d3 from "d3";
import { Colors } from "@blueprintjs/core";

const Lasso = () => {
  const dispatch = d3.dispatch("start", "end", "cancel");

  const lasso = (svg) => {
    let lassoPolygon;
    let lassoPath;
    let closePath;
    let lassoInProgress;

    const polygonToPath = (polygon) =>
      `M${polygon.map((d) => d.join(",")).join("L")}`;

    const distance = (pt1, pt2) =>
      Math.sqrt((pt2[0] - pt1[0]) ** 2 + (pt2[1] - pt1[1]) ** 2);

    // distance last point has to be to first point before it auto closes when mouse is released
    const closeDistance = 75;
    const lassoPathColor = Colors.BLUE5;

    const handleDragStart = () => {
      lassoPolygon = [d3.mouse(svg.node())]; // current x y of mouse within element

      if (lassoPath) {
        // If the existing path is in progress
        if (lassoInProgress) {
          // cancel the existing lasso
          handleCancel();
          // Don't continue with current drag start
          return;
        }

        lassoPath.remove();
      }
      // We're starting a new drag
      lassoInProgress = true;

      lassoPath = g
        .append("path")
        .attr("data-testid", "lasso-element")
        .attr("fill-opacity", 0.1)
        .attr("stroke-dasharray", "3, 3");

      closePath = g
        .append("line")
        .attr("x2", lassoPolygon[0][0])
        .attr("y2", lassoPolygon[0][1])
        .attr("stroke-dasharray", "3, 3");

      dispatch.call("start", lasso, lassoPolygon);
    };

    const handleDrag = () => {
      const point = d3.mouse(svg.node());
      lassoPolygon.push(point);
      lassoPath.attr("d", polygonToPath(lassoPolygon));

      // indicate if we are within closing distance
      if (
        distance(lassoPolygon[0], lassoPolygon[lassoPolygon.length - 1]) <
        closeDistance
      ) {
        const closePathColor = Colors.GREEN5;
        closePath
          .attr("x1", point[0])
          .attr("y1", point[1])
          .attr("opacity", 1)
          .attr("stroke", closePathColor)
          .attr("fill", closePathColor);
        lassoPath.attr("stroke", closePathColor).attr("fill", closePathColor);
      } else {
        closePath.attr("opacity", 0);
        lassoPath.attr("stroke", lassoPathColor).attr("fill", lassoPathColor);
      }
    };

    const handleCancel = () => {
      lassoPath.remove();
      closePath = closePath?.remove();
      lassoPath = null;
      lassoPolygon = null;
      closePath = null;
      dispatch.call("cancel");
    };

    const handleDragEnd = () => {
      // remove the close path
      closePath.remove();
      closePath = null;

      // successfully closed
      if (
        distance(lassoPolygon[0], lassoPolygon[lassoPolygon.length - 1]) <
        closeDistance
      ) {
        lassoInProgress = false;

        lassoPath.attr("d", `${polygonToPath(lassoPolygon)}Z`);
        dispatch.call("end", lasso, lassoPolygon);

        // otherwise cancel
      } else {
        handleCancel();
      }
    };

    // append a <g> with a rect
    const g = svg.append("g").attr("class", "lasso-group");
    const bbox = svg.node().getBoundingClientRect();
    const area = g
      .append("rect")
      .attr("width", bbox.width)
      .attr("height", bbox.height)
      .attr("opacity", 0);

    const drag = d3
      .drag()
      .on("start", handleDragStart)
      .on("drag", handleDrag)
      .on("end", handleDragEnd);

    area.call(drag);

    lasso.reset = () => {
      if (lassoPath) {
        lassoPath.remove();
        lassoPath = null;
      }

      lassoPolygon = null;
      if (closePath) {
        closePath.remove();
        closePath = null;
      }
    };

    lasso.move = (polygon) => {
      if (polygon !== lassoPolygon || polygon.length !== lassoPolygon.length) {
        lasso.reset();

        lassoPolygon = polygon;
        lassoPath = g
          .append("path")
          .attr("data-testid", "lasso-element")
          .attr("fill", lassoPathColor)
          .attr("fill-opacity", 0.1)
          .attr("stroke", lassoPathColor)
          .attr("stroke-dasharray", "3, 3");

        lassoPath.attr("d", `${polygonToPath(lassoPolygon)}Z`);
      }
    };
  };

  lasso.on = (type, callback) => {
    dispatch.on(type, callback);
    return lasso;
  };

  return lasso;
};

export default Lasso;
