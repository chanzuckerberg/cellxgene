import React, { useEffect, useRef } from "react";
import { connect } from "react-redux";
import { Button, ButtonGroup } from "@blueprintjs/core";
import _regl from "regl";
import * as d3 from "d3";
import { mat3 } from "gl-matrix";
import memoize from "memoize-one";
import Async from "react-async";

import * as globals from "../../globals";
// import setupScatterplot from "./setupScatterplot";
import styles from "./scatterplot.css";
import _drawPoints from "./drawPointsRegl";
import { margin, width, height } from "./util";
import {
  createColorTable,
  createColorQuery,
} from "../../util/stateManager/colorHelpers";
import renderThrottle from "../../util/renderThrottle";
import shallowEqual from "../../util/shallowEqual";

const flagSelected = 1;
const flagNaN = 2;
const flagHighlight = 4;

function createProjectionTF(viewportWidth, viewportHeight) {
  /*
  the projection transform accounts for the screen size & other layout
  */
  const m = mat3.create();
  return mat3.projection(m, viewportWidth, viewportHeight);
}

@connect((state) => {
  const { obsCrossfilter: crossfilter } = state;
  const { scatterplotXXaccessor, scatterplotYYaccessor } = state.controls;

  return {
    annoMatrix: state.annoMatrix,
    colors: state.colors,
    pointDilation: state.pointDilation,

    // Accessors are var/gene names (strings)
    scatterplotXXaccessor,
    scatterplotYYaccessor,

    differential: state.differential,
    crossfilter,
  };
})
class Scatterplot extends React.PureComponent {
  static createReglState(canvas) {
    /*
    Must be created for each canvas
    */
    // setup canvas, webgl draw function and camera
    const regl = _regl(canvas);
    const drawPoints = _drawPoints(regl);

    // preallocate webgl buffers
    const pointBuffer = regl.buffer();
    const colorBuffer = regl.buffer();
    const flagBuffer = regl.buffer();

    return {
      regl,
      drawPoints,
      pointBuffer,
      colorBuffer,
      flagBuffer,
    };
  }

  static watchAsync(props, prevProps) {
    return !shallowEqual(props.watchProps, prevProps.watchProps);
  }

  static getScale(col, rangeMin, rangeMax) {
    if (!col) return null;
    const { min, max } = col.summarize();
    return d3.scaleLinear().domain([min, max]).range([rangeMin, rangeMax]);
  }

  computePointPositions = memoize((X, Y, xScale, yScale) => {
    const positions = new Float32Array(2 * X.length);
    for (let i = 0, len = X.length; i < len; i += 1) {
      positions[2 * i] = xScale(X[i]);
      positions[2 * i + 1] = yScale(Y[i]);
    }
    return positions;
  });

  computePointColors = memoize((rgb) => {
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
    (crossfilter, _flagSelected, _flagUnselected) => {
      const x = crossfilter.fillByIsSelected(
        new Float32Array(crossfilter.size()),
        _flagSelected,
        _flagUnselected
      );
      return x;
    }
  );

  computePointFlags = memoize(
    (crossfilter, colorByData, pointDilationData, pointDilationLabel) => {
      /*
      We communicate with the shader using three flags:
      - isNaN -- the value is a NaN. Only makes sense when we have a colorAccessor
      - isSelected -- the value is selected
      - isHightlighted -- the value is highlighted in the UI (orthogonal from selection highlighting)

      Due to constraints in webgl vertex shader attributes, these are encoded in a float, "kinda"
      like bitmasks.

      We also have separate code paths for generating flags for categorical and
      continuous metadata, as they rely on different tests, and some of the flags
      (eg, isNaN) are meaningless in the face of categorical metadata.
      */

      const flags = this.computeSelectedFlags(
        crossfilter,
        flagSelected,
        0
      ).slice();

      if (colorByData || pointDilationData) {
        for (let i = 0, len = flags.length; i < len; i += 1) {
          if (pointDilationData) {
            flags[i] +=
              pointDilationData[i] === pointDilationLabel ? flagHighlight : 0;
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
    this.axes = false;
    this.reglCanvas = null;
    this.renderCache = null;
    this.state = {
      regl: null,
      // svg: null,
      drawPoints: null,
      minimized: null,
      viewport: {
        height: null,
        width: null,
      },
      projectionTF: null,
    };
  }

  componentDidMount() {
    // Create render transform
    const projectionTF = createProjectionTF(
      this.reglCanvas.width,
      this.reglCanvas.height
    );

    window.addEventListener("resize", this.handleResize);
    const viewport = this.getViewportDimensions();

    this.setState({
      projectionTF,
      viewport,
    });
  }

  componentWillUnmount() {
    window.removeEventListener("resize", this.updateViewportDimensions);
  }

  setReglCanvas = (canvas) => {
    this.reglCanvas = canvas;
    this.setState({
      ...Scatterplot.createReglState(canvas),
    });
  };

  getViewportDimensions = () => {
    return {
      viewport: {
        height: window.height,
        width: window.width,
      },
    };
  };

  handleResize = () => {
    const { state } = this.state;
    const viewport = this.getViewportDimensions();
    this.setState({
      ...state,
      viewport,
    });
  };

  updateViewportDimensions = () => {
    this.setState(this.getViewportDimensions());
  };

  fetchAsyncProps = async (props) => {
    const {
      scatterplotXXaccessor,
      scatterplotYYaccessor,
      colors: colorsProp,
      crossfilter,
      pointDilation,
    } = props.watchProps;

    const [
      expressionXDf,
      expressionYDf,
      colorDf,
      pointDilationDf,
    ] = await this.fetchData(
      scatterplotXXaccessor,
      scatterplotYYaccessor,
      colorsProp,
      pointDilation
    );
    const colorTable = this.updateColorTable(colorsProp, colorDf);

    const xCol = expressionXDf.icol(0);
    const yCol = expressionYDf.icol(0);
    const xScale = Scatterplot.getScale(xCol, 0, width);
    const yScale = Scatterplot.getScale(yCol, height, 0);
    const positions = this.computePointPositions(
      xCol.asArray(),
      yCol.asArray(),
      xScale,
      yScale
    );

    const colors = this.computePointColors(colorTable.rgb);

    const { colorAccessor } = colorsProp;
    const colorByData = colorDf?.col(colorAccessor)?.asArray();
    const {
      metadataField: pointDilationCategory,
      categoryField: pointDilationLabel,
    } = pointDilation;
    const pointDilationData = pointDilationDf
      ?.col(pointDilationCategory)
      ?.asArray();
    const flags = this.computePointFlags(
      crossfilter,
      colorByData,
      pointDilationData,
      pointDilationLabel
    );

    return {
      positions,
      colors,
      flags,
      width,
      height,
      xScale,
      yScale,
    };
  };

  createXQuery(geneName) {
    const { annoMatrix } = this.props;
    const { schema } = annoMatrix;
    const varIndex = schema?.annotations?.var?.index;
    if (!varIndex) return null;
    return [
      "X",
      {
        field: "var",
        column: varIndex,
        value: geneName,
      },
    ];
  }

  createColorByQuery(colors) {
    const { annoMatrix } = this.props;
    const { schema } = annoMatrix;
    const { colorMode, colorAccessor } = colors;
    return createColorQuery(colorMode, colorAccessor, schema);
  }

  updateColorTable(colors, colorDf) {
    /* update color table state */
    const { annoMatrix } = this.props;
    const { schema } = annoMatrix;
    const { colorAccessor, userColors, colorMode } = colors;
    return createColorTable(
      colorMode,
      colorAccessor,
      colorDf,
      schema,
      userColors
    );
  }

  async fetchData(
    scatterplotXXaccessor,
    scatterplotYYaccessor,
    colors,
    pointDilation
  ) {
    const { annoMatrix } = this.props;
    const { metadataField: pointDilationAccessor } = pointDilation;

    const promises = [];
    // X and Y dimensions
    promises.push(
      annoMatrix.fetch(...this.createXQuery(scatterplotXXaccessor))
    );
    promises.push(
      annoMatrix.fetch(...this.createXQuery(scatterplotYYaccessor))
    );

    // color
    const query = this.createColorByQuery(colors);
    if (query) {
      promises.push(annoMatrix.fetch(...query));
    } else {
      promises.push(Promise.resolve(null));
    }

    // point highlighting
    if (pointDilationAccessor) {
      promises.push(annoMatrix.fetch("obs", pointDilationAccessor));
    } else {
      promises.push(Promise.resolve(null));
    }

    return Promise.all(promises);
  }

  renderCanvas = renderThrottle(() => {
    const {
      regl,
      drawPoints,
      colorBuffer,
      pointBuffer,
      flagBuffer,
      projectionTF,
    } = this.state;
    this.renderPoints(
      regl,
      drawPoints,
      flagBuffer,
      colorBuffer,
      pointBuffer,
      projectionTF
    );
  });

  updateReglAndRender(newRenderCache) {
    const { positions, colors, flags } = newRenderCache;
    this.renderCache = newRenderCache;
    const { pointBuffer, colorBuffer, flagBuffer } = this.state;
    pointBuffer({ data: positions, dimension: 2 });
    colorBuffer({ data: colors, dimension: 3 });
    flagBuffer({ data: flags, dimension: 1 });
    this.renderCanvas();
  }

  renderPoints(
    regl,
    drawPoints,
    flagBuffer,
    colorBuffer,
    pointBuffer,
    projectionTF
  ) {
    const { annoMatrix } = this.props;
    if (!this.reglCanvas || !annoMatrix) return;

    const { schema } = annoMatrix;
    const { viewport } = this.state;
    regl.poll();
    regl.clear({
      depth: 1,
      color: [1, 1, 1, 1],
    });
    drawPoints({
      flag: flagBuffer,
      color: colorBuffer,
      position: pointBuffer,
      projection: projectionTF,
      count: annoMatrix.nObs,
      nPoints: schema.dataframe.nObs,
      minViewportDimension: Math.min(
        viewport.width - globals.leftSidebarWidth || width,
        viewport.height || height
      ),
    });
    regl._gl.flush();
  }

  render() {
    const {
      dispatch,
      annoMatrix,
      scatterplotXXaccessor,
      scatterplotYYaccessor,
      colors,
      crossfilter,
      pointDilation,
    } = this.props;
    const { minimized, status, regl, viewport } = this.state;

    if (status === "error") return null;
    if (regl) {
      this.renderCanvas();
    }

    return (
      <div
        style={{
          position: "fixed",
          bottom: minimized ? -height + -margin.top - 2 : 0,
          borderRadius: "3px 3px 0px 0px",
          left: globals.leftSidebarWidth + globals.scatterplotMarginLeft,
          padding: "0px 20px 20px 0px",
          background: "white",
          /* x y blur spread color */
          boxShadow: "0px 0px 6px 2px rgba(153,153,153,0.4)",
          zIndex: 2,
        }}
        id="scatterplot_wrapper"
      >
        <ButtonGroup
          style={{
            position: "absolute",
            right: 5,
            top: 5,
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
                type: "clear scatterplot",
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
            height: `${height + margin.top + margin.bottom}px`,
          }}
        >
          <canvas
            width={width}
            height={height}
            data-testid="scatterplot"
            style={{
              marginLeft: margin.left,
              marginTop: margin.top,
            }}
            ref={this.setReglCanvas}
          />
          <Async
            watchFn={Scatterplot.watchAsync}
            promiseFn={this.fetchAsyncProps}
            watchProps={{
              annoMatrix,
              scatterplotXXaccessor,
              scatterplotYYaccessor,
              colors,
              crossfilter,
              pointDilation,
              viewport,
            }}
          >
            <Async.Pending>Loading...</Async.Pending>
            <Async.Rejected>{(error) => error.message}</Async.Rejected>
            <Async.Fulfilled>
              {(asyncProps) => {
                if (regl && !shallowEqual(asyncProps, this.renderCache)) {
                  this.updateReglAndRender(asyncProps);
                }
                return (
                  <ScatterplotAxis
                    width={width}
                    height={height}
                    margin={margin}
                    scatterplotYYaccessor={scatterplotXXaccessor}
                    scatterplotXXaccessor={scatterplotYYaccessor}
                    xScale={asyncProps.xScale}
                    yScale={asyncProps.yScale}
                  />
                );
              }}
            </Async.Fulfilled>
          </Async>
        </div>
      </div>
    );
  }
}

export default Scatterplot;

const ScatterplotAxis = React.memo(
  ({ scatterplotYYaccessor, scatterplotXXaccessor, xScale, yScale }) => {
    /*
    Axis for the scatterplot, rendered with SVG/D3.  Props:
      * scatterplotXXaccessor - name of X axis
      * scatterplotXXaccessor - name of Y axis
      * xScale - D3 scale for X axis (domain to range)
      * yScale - D3 scale for Y axis (domain to range)

    This also relies on the GLOBAL width/height/margin constants.  If those become
    become variables, may need to add the params.
    */

    const svgRef = useRef(null);

    useEffect(() => {
      if (!svgRef.current) return;
      const svg = d3.select(svgRef.current);

      svg.selectAll("*").remove();

      // the axes are much cleaner and easier now. No need to rotate and orient
      // the axis, just call axisBottom, axisLeft etc.
      const xAxis = d3.axisBottom().ticks(7).scale(xScale);
      const yAxis = d3.axisLeft().ticks(7).scale(yScale);

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
    }, [scatterplotXXaccessor, scatterplotYYaccessor, xScale, yScale]);

    return (
      <svg
        width={width + margin.left + margin.right}
        height={height + margin.top + margin.bottom}
        data-testid="scatterplot-svg"
      >
        <g ref={svgRef} transform={`translate(${margin.left},${margin.top})`} />
      </svg>
    );
  }
);
