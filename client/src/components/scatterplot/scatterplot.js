import React from "react";
import { connect } from "react-redux";
import { Button, ButtonGroup } from "@blueprintjs/core";
import _regl from "regl";
import * as d3 from "d3";
import { mat3 } from "gl-matrix";
import memoize from "memoize-one";
import { isTypedArray } from "../../util/typeHelpers";

import * as globals from "../../globals";
import setupScatterplot from "./setupScatterplot";
import styles from "./scatterplot.css";
import _drawPoints from "./drawPointsRegl";
import { margin, width, height } from "./util";
import finiteExtent from "../../util/finiteExtent";

function createProjectionTF(viewportWidth, viewportHeight) {
  /*
  the projection transform accounts for the screen size & other layout
  */
  const m = mat3.create();
  return mat3.projection(m, viewportWidth, viewportHeight);
}

@connect(state => {
  const { world, crossfilter, universe } = state;
  const { scatterplotXXaccessor, scatterplotYYaccessor } = state.controls;
  const expressionX = scatterplotXXaccessor
    ? world.varData.col(scatterplotXXaccessor)?.asArray()
    : null;
  const expressionY = scatterplotYYaccessor
    ? world.varData.col(scatterplotYYaccessor)?.asArray()
    : null;

  return {
    world,
    universe,

    colorRGB: state.colors.rgb,
    colorScale: state.colors.scale,
    colorAccessor: state.colors.colorAccessor,

    centroidLabel: state.centroidLabel,

    // Accessors are var/gene names (strings)
    scatterplotXXaccessor,
    scatterplotYYaccessor,
    opacityForDeselectedCells: state.controls.opacityForDeselectedCells,

    differential: state.differential,

    expressionX,
    expressionY,

    crossfilter,

    responsive: state.responsive
  };
})
class Scatterplot extends React.PureComponent {
  computePointPositions = memoize((X, Y, xScale, yScale) => {
    const positions = new Float32Array(2 * X.length);
    for (let i = 0, len = X.length; i < len; i += 1) {
      positions[2 * i] = xScale(X[i]);
      positions[2 * i + 1] = yScale(Y[i]);
    }
    return positions;
  });

  computePointColors = memoize(rgb => {
    /*
    compute webgl colors for each point
    */
    const colors = new Float32Array(3 * rgb.length);
    for (let i = 0, len = rgb.length; i < len; i += 1) {
      colors.set(rgb[i], 3 * i);
    }
    return colors;
  });

  computeSelectedFlags = memoize(
    (crossfilter, flagSelected, flagUnselected) => {
      const x = crossfilter.fillByIsSelected(
        new Float32Array(crossfilter.size()),
        flagSelected,
        flagUnselected
      );
      return x;
    }
  );

  computePointFlags = memoize(
    (world, crossfilter, colorAccessor, centroidLabel) => {
      const flagSelected = 1;
      const flagNaN = 2;
      const flagHighlight = 4;

      const flags = this.computeSelectedFlags(
        crossfilter,
        flagSelected,
        0
      ).slice();

      const { metadataField, categoryField } = centroidLabel;
      const highlightData = metadataField
        ? world.obsAnnotations.col(metadataField)?.asArray()
        : null;
      const colorByColumn = colorAccessor
        ? world.obsAnnotations.col(colorAccessor)?.asArray() ||
          world.varData.col(colorAccessor)?.asArray()
        : null;
      const colorByData =
        colorByColumn && isTypedArray(colorByColumn) ? colorByColumn : null;

      if (colorByData || highlightData) {
        for (let i = 0, len = flags.length; i < len; i += 1) {
          if (highlightData) {
            flags[i] += highlightData[i] === categoryField ? flagHighlight : 0;
          }
          if (colorByData) {
            flags[i] += Number.isFinite(colorByData[i]) ? 0 : flagNaN;
          }
        }
      }
      return flags;
    }
  );

  constructor(props) {
    super(props);
    this.count = 0;
    this.axes = false;
    this.renderCache = {
      positions: null,
      colors: null,
      flags: null,
      xScale: null,
      yScale: null
    };
    this.state = {
      svg: null,
      minimized: null
    };
  }

  componentDidMount() {
    const { svg } = setupScatterplot(width, height, margin);
    let scales;
    const { expressionX, expressionY } = this.props;

    if (svg && expressionX && expressionY) {
      scales = Scatterplot.setupScales(expressionX, expressionY);
      this.drawAxesSVG(scales.xScale, scales.yScale, svg);
      this.renderCache = { ...this.renderCache, ...scales };
    }

    const regl = _regl(this.reglCanvas);
    const drawPoints = _drawPoints(regl);

    // Create render transform
    const projectionTF = createProjectionTF(
      this.reglCanvas.width,
      this.reglCanvas.height
    );

    // preallocate buffers
    const pointBuffer = regl.buffer();
    const colorBuffer = regl.buffer();
    const flagBuffer = regl.buffer();

    this.renderPoints(
      regl,
      drawPoints,
      flagBuffer,
      colorBuffer,
      pointBuffer,
      projectionTF
    );

    this.setState({
      regl,
      flagBuffer,
      pointBuffer,
      colorBuffer,
      svg,
      drawPoints,
      projectionTF
    });
  }

  componentDidUpdate(prevProps) {
    const {
      world,
      crossfilter,
      scatterplotXXaccessor,
      scatterplotYYaccessor,
      expressionX,
      expressionY,
      colorRGB,
      colorAccessor,
      centroidLabel
    } = this.props;
    const {
      regl,
      pointBuffer,
      colorBuffer,
      flagBuffer,
      svg,
      drawPoints,
      projectionTF
    } = this.state;

    if (
      scatterplotXXaccessor !== prevProps.scatterplotXXaccessor ||
      scatterplotYYaccessor !== prevProps.scatterplotYYaccessor ||
      world !== prevProps.world // shape or clip of world changed
    ) {
      const scales = Scatterplot.setupScales(expressionX, expressionY);
      this.drawAxesSVG(scales.xScale, scales.yScale, svg);
      this.renderCache = { ...this.renderCache, ...scales };
    }

    if (world && regl) {
      const { renderCache } = this;
      const { xScale, yScale } = this.renderCache;
      let needsRepaint = false;

      const newPositions = this.computePointPositions(
        expressionX,
        expressionY,
        xScale,
        yScale
      );
      if (renderCache.positions !== newPositions) {
        renderCache.positions = newPositions;
        pointBuffer({ data: renderCache.positions, dimension: 2 });
        needsRepaint = true;
      }

      /* colors for each point */
      const newColors = this.computePointColors(colorRGB);
      if (renderCache.colors !== newColors) {
        renderCache.colors = newColors;
        colorBuffer({ data: renderCache.colors, dimension: 3 });
        needsRepaint = true;
      }

      const newFlags = this.computePointFlags(
        world,
        crossfilter,
        colorAccessor,
        centroidLabel
      );
      if (renderCache.flags !== newFlags) {
        renderCache.flags = newFlags;
        flagBuffer({ data: renderCache.flags, dimension: 1 });
        needsRepaint = true;
      }

      this.count = expressionX.length;

      if (needsRepaint) {
        this.renderPoints(
          regl,
          drawPoints,
          flagBuffer,
          colorBuffer,
          pointBuffer,
          projectionTF
        );
      }
    }
  }

  static setupScales(expressionX, expressionY) {
    const xScale = d3
      .scaleLinear()
      .domain(finiteExtent(expressionX))
      .range([0, width]);
    const yScale = d3
      .scaleLinear()
      .domain(finiteExtent(expressionY))
      .range([height, 0]);

    return {
      xScale,
      yScale
    };
  }

  drawAxesSVG(xScale, yScale, svg) {
    const { scatterplotYYaccessor, scatterplotXXaccessor } = this.props;
    svg.selectAll("*").remove();

    // the axes are much cleaner and easier now. No need to rotate and orient
    // the axis, just call axisBottom, axisLeft etc.
    const xAxis = d3
      .axisBottom()
      .ticks(7)
      .scale(xScale);

    const yAxis = d3
      .axisLeft()
      .ticks(7)
      .scale(yScale);

    // adding axes is also simpler now, just translate x-axis to (0,height)
    // and it's alread defined to be a bottom axis.
    svg
      .append("g")
      .attr("transform", `translate(0,${height})`)
      .attr("class", "x axis")
      .call(xAxis);

    // y-axis is translated to (0,0)
    svg
      .append("g")
      .attr("transform", "translate(0,0)")
      .attr("class", "y axis")
      .call(yAxis);

    // adding label. For x-axis, it's at (10, 10), and for y-axis at (width, height-10).
    svg
      .append("text")
      .attr("x", 10)
      .attr("y", 10)
      .attr("class", "label")
      .style("font-style", "italic")
      .text(scatterplotYYaccessor);

    svg
      .append("text")
      .attr("x", width)
      .attr("y", height - 10)
      .attr("text-anchor", "end")
      .attr("class", "label")
      .style("font-style", "italic")
      .text(scatterplotXXaccessor);
  }

  renderPoints(
    regl,
    drawPoints,
    flagBuffer,
    colorBuffer,
    pointBuffer,
    projectionTF
  ) {
    if (!this.reglCanvas) return;
    const { universe, responsive } = this.props;
    // The viewport dimension is used to scale points, so we want to pass
    // the dimension of the MAIN viewport, not the scatterplot viewport.
    // Slightly hacky, but we want all points to scale uniformly.  Perhaps
    // this should move to the redux state and be shared?
    const { width: cvWidth, height: cvHeight } = responsive;
    regl.poll();
    regl.clear({
      depth: 1,
      color: [1, 1, 1, 1]
    });
    drawPoints({
      flag: flagBuffer,
      color: colorBuffer,
      position: pointBuffer,
      projection: projectionTF,
      count: this.count,
      nPoints: universe.nObs,
      minViewportDimension: Math.min(
        cvWidth - globals.leftSidebarWidth || width,
        cvHeight || height
      )
    });
    regl._gl.flush();
  }

  render() {
    const { dispatch } = this.props;
    const { minimized } = this.state;

    return (
      <div
        style={{
          position: "fixed",
          bottom: minimized ? -height + -margin.top : 0,
          borderRadius: "3px 3px 0px 0px",
          left: globals.leftSidebarWidth + globals.scatterplotMarginLeft,
          padding: "0px 20px 20px 0px",
          backgroundColor: "white",
          /* x y blur spread color */
          boxShadow: "0px 0px 6px 2px rgba(153,153,153,0.4)"
        }}
        id="scatterplot_wrapper"
      >
        <ButtonGroup
          style={{
            position: "absolute",
            right: 5,
            top: 5
          }}
        >
          <Button
            type="button"
            minimal
            onClick={() => {
              this.setState({ minimized: !minimized });
            }}
          >
            {minimized ? "show scatterplot" : "hide"}
          </Button>
          <Button
            type="button"
            minimal
            data-testid="clear-scatterplot"
            onClick={() =>
              dispatch({
                type: "clear scatterplot"
              })
            }
          >
            remove
          </Button>
        </ButtonGroup>
        <div
          className={styles.scatterplot}
          id="scatterplot"
          style={{
            width: `${width + margin.left + margin.right}px`,
            height: `${height + margin.top + margin.bottom}px`
          }}
        >
          <canvas
            width={width}
            height={height}
            data-testid="scatterplot"
            style={{
              marginLeft: margin.left,
              marginTop: margin.top
            }}
            ref={canvas => {
              this.reglCanvas = canvas;
            }}
          />
        </div>
      </div>
    );
  }
}

export default Scatterplot;
