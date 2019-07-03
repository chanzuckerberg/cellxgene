// https://bl.ocks.org/pbeshai/8008075f9ce771ee8be39e8c38907570

import * as d3 from "d3";

const Lasso = () => {
  const dispatch = d3.dispatch("start", "end", "cancel");

  const polygonToPath = polygon =>
    `M${polygon.map(d => d.join(",")).join("L")}`;

  const distance = (pt1, pt2) =>
    Math.sqrt((pt2[0] - pt1[0]) ** 2 + (pt2[1] - pt1[1]) ** 2);

  // distance last point has to be to first point before it auto closes when mouse is released
  const closeDistance = 75;

  const lasso = svg => {
    let lassoPolygon;
    let lassoPath;
    let closePath;

    const handleDragStart = () => {
      lassoPolygon = [d3.mouse(svg.node())]; // current x y of mouse within element

      if (lassoPath) {
        lassoPath.remove();
      }

      lassoPath = g
        .append("path")
        .attr("data-testid", "lasso-element")
        .attr("fill", "#0bb")
        .attr("fill-opacity", 0.1)
        .attr("stroke", "#0bb")
        .attr("stroke-dasharray", "3, 3");

      closePath = g
        .append("line")
        .attr("x2", lassoPolygon[0][0])
        .attr("y2", lassoPolygon[0][1])
        .attr("stroke", "#0bb")
        .attr("stroke-dasharray", "3, 3")
        .attr("opacity", 0);

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
        closePath
          .attr("x1", point[0])
          .attr("y1", point[1])
          .attr("opacity", 1);
      } else {
        closePath.attr("opacity", 0);
      }
    };

    const handleDragEnd = () => {
      // remove the close path
      closePath.remove();
      closePath = null;

      // succesfully closed
      if (
        distance(lassoPolygon[0], lassoPolygon[lassoPolygon.length - 1]) <
        closeDistance
      ) {
        lassoPath.attr("d", `${polygonToPath(lassoPolygon)}Z`);
        dispatch.call("end", lasso, lassoPolygon);

        // otherwise cancel
      } else {
        lassoPath.remove();
        lassoPath = null;
        lassoPolygon = null;
        dispatch.call("cancel");
      }
    };

    // append a <g> with a rect
    const g = svg.append("g").attr("class", "lasso-group");
    const bbox = svg.node().getBoundingClientRect();
    const area = g
      .append("rect")
      .attr("width", bbox.width)
      .attr("height", bbox.height)
      .attr("fill", "tomato")
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

    lasso.move = polygon => {
      if (polygon !== lassoPolygon || polygon.length !== lassoPolygon.length) {
        lasso.reset();

        lassoPolygon = polygon;
        lassoPath = g
          .append("path")
          .attr("data-testid", "lasso-element")
          .attr("fill", "#0bb")
          .attr("fill-opacity", 0.1)
          .attr("stroke", "#0bb")
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
