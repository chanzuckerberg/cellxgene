// jshint esversion: 6
// https://bl.ocks.org/Jverma/076377dd0125b1a508621441752735fc
// https://peterbeshai.com/scatterplot-in-d3-with-voronoi-interaction.html

import React from "react";
import _ from "lodash";
import { connect } from "react-redux";
import { Button, ButtonGroup } from "@blueprintjs/core";
import _regl from "regl";
import * as d3 from "d3";
import * as globals from "../../globals";

import _camera from "../../util/camera";

import setupScatterplot from "./setupScatterplot";
import styles from "./scatterplot.css";

import _drawPoints from "./drawPointsRegl";
import scaleLinear from "../../util/scaleLinear";

import { margin, width, height } from "./util";
import finiteExtent from "../../util/finiteExtent";

@connect(state => {
  const { world, crossfilter } = state;
  const { scatterplotXXaccessor, scatterplotYYaccessor } = state.controls;
  const expressionX =
    world &&
    scatterplotXXaccessor &&
    world.varData.hasCol(scatterplotXXaccessor)
      ? world.varData.col(scatterplotXXaccessor).asArray()
      : null;
  const expressionY =
    world &&
    scatterplotYYaccessor &&
    world.varData.hasCol(scatterplotYYaccessor)
      ? world.varData.col(scatterplotYYaccessor).asArray()
      : null;

  return {
    world,

    colorRGB: state.colors.rgb,
    colorScale: state.colors.scale,
    colorAccessor: state.colors.colorAccessor,

    // Accessors are var/gene names (strings)
    scatterplotXXaccessor,
    scatterplotYYaccessor,
    opacityForDeselectedCells: state.controls.opacityForDeselectedCells,

    differential: state.differential,

    expressionX,
    expressionY,

    crossfilter
  };
})
class Scatterplot extends React.Component {
  constructor(props) {
    super(props);
    this.count = 0;
    this.axes = false;
    this.renderCache = {
      positions: null,
      colors: null,
      sizes: null,
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

    const camera = _camera(this.reglCanvas, { scale: true, rotate: false });
    const regl = _regl(this.reglCanvas);

    const drawPoints = _drawPoints(regl);

    // preallocate buffers
    const pointBuffer = regl.buffer();
    const colorBuffer = regl.buffer();
    const sizeBuffer = regl.buffer();

    const reglRender = regl.frame(() => {
      this.reglDraw(
        regl,
        drawPoints,
        sizeBuffer,
        colorBuffer,
        pointBuffer,
        camera
      );
      camera.tick();
    });

    this.reglRenderState = "rendering";

    this.setState({
      regl,
      sizeBuffer,
      pointBuffer,
      colorBuffer,
      svg,
      reglRender,
      camera,
      drawPoints
    });
  }

  componentDidUpdate(prevProps) {
    console.log("scatterplot::componentDidUpdate");
    const {
      world,
      crossfilter,
      scatterplotXXaccessor,
      scatterplotYYaccessor,
      expressionX,
      expressionY,
      colorRGB
    } = this.props;
    const {
      reglRender,
      regl,
      pointBuffer,
      colorBuffer,
      sizeBuffer,
      svg,
      drawPoints,
      camera
    } = this.state;

    if (
      scatterplotXXaccessor !== prevProps.scatterplotXXaccessor || // was CLU now FTH1 etc
      scatterplotYYaccessor !== prevProps.scatterplotYYaccessor || // was CLU now FTH1 etc
      world !== prevProps.world // shape or clip of world changed
    ) {
      const scales = Scatterplot.setupScales(expressionX, expressionY);
      this.drawAxesSVG(scales.xScale, scales.yScale, svg);
      this.renderCache = { ...this.renderCache, ...scales };
    }

    if (reglRender && this.reglRenderState === "rendering") {
      reglRender.cancel();
      this.reglRenderState = "paused";
    }

    if (
      world &&
      regl &&
      pointBuffer &&
      colorBuffer &&
      sizeBuffer &&
      expressionX &&
      expressionY &&
      scatterplotXXaccessor &&
      scatterplotYYaccessor
    ) {
      const { renderCache } = this;
      const { xScale, yScale } = this.renderCache;
      const cellCount = expressionX.length;

      // Points change when expressionX or expressionY change.
      if (
        !renderCache.positions ||
        expressionX !== prevProps.expressionX ||
        expressionY !== prevProps.expressionY
      ) {
        if (!renderCache.positions) {
          renderCache.positions = new Float32Array(2 * cellCount);
        }
        const glScaleX = scaleLinear([0, width], [-0.95, 0.95]);
        const glScaleY = scaleLinear([0, height], [-1, 1]);
        for (let i = 0, { positions } = renderCache; i < cellCount; i += 1) {
          positions[2 * i] = glScaleX(xScale(expressionX[i]));
          positions[2 * i + 1] = glScaleY(yScale(expressionY[i]));
        }
        pointBuffer({ data: renderCache.positions, dimension: 2 });
      }

      // Colors for each point - change only when props.colorsRGB change.
      if (!renderCache.colors || colorRGB !== prevProps.colorRGB) {
        if (!renderCache.colors) {
          renderCache.colors = new Float32Array(3 * cellCount);
        }
        for (let i = 0, { colors } = renderCache; i < cellCount; i += 1) {
          colors.set(colorRGB[i], 3 * i);
        }
        colorBuffer({ data: renderCache.colors, dimension: 3 });
      }

      // Sizes for each point - updates are triggered only when selected
      // obs change
      if (!renderCache.sizes || crossfilter !== prevProps.crossfilter) {
        if (!renderCache.sizes) {
          renderCache.sizes = new Float32Array(cellCount);
        }
        crossfilter.fillByIsSelected(renderCache.sizes, 4, 0.2);
        sizeBuffer({ data: renderCache.sizes, dimension: 1 });
      }

      this.count = cellCount;

      regl._refresh();
      this.reglDraw(
        regl,
        drawPoints,
        sizeBuffer,
        colorBuffer,
        pointBuffer,
        camera
      );
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

  reglDraw(regl, drawPoints, sizeBuffer, colorBuffer, pointBuffer, camera) {
    regl.clear({
      depth: 1,
      color: [1, 1, 1, 1]
    });

    drawPoints({
      size: sizeBuffer,
      distance: camera.distance,
      color: colorBuffer,
      position: pointBuffer,
      count: this.count,
      view: camera.view()
    });
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
              marginLeft: margin.left - 7,
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
