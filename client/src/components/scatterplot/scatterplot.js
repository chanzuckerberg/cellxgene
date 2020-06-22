import React from "react";
import { connect } from "react-redux";
import { Button, ButtonGroup } from "@blueprintjs/core";
import _regl from "regl";
import * as d3 from "d3";
import { mat3 } from "gl-matrix";
import memoize from "memoize-one";

import * as globals from "../../globals";
import setupScatterplot from "./setupScatterplot";
import styles from "./scatterplot.css";
import _drawPoints from "./drawPointsRegl";
import { margin, width, height } from "./util";
import {
  createColorTable,
  createColorQuery,
} from "../../util/stateManager/colorHelpers";
import renderThrottle from "../../util/renderThrottle";

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
  static drawAxesSVG(
    xScale,
    yScale,
    svg,
    scatterplotYYaccessor,
    scatterplotXXaccessor
  ) {
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
    this.renderCache = {
      // cached state that doesn't trigger an update
      positions: null,
      colors: null,
      flags: null,
      xScale: null,
      yScale: null,
    };
    this.state = {
      regl: null,
      svg: null,
      drawPoints: null,
      minimized: null,
      viewport: {
        height: null,
        width: null,
      },
      projectionTF: null,

      // component rendering derived state - these must stay synchronized
      // with the reducer state they were generated from.
      scatterplotState: {
        scatterplotXXaccessor: null,
        scatterplotYYaccessor: null,
        expressionXDf: null,
        expressionYDf: null,
      },
      colorState: {
        colors: null,
        colorDf: null,
        colorTable: null,
      },
      pointDilationState: {
        pointDilation: null,
        pointDilationDf: null,
      },
    };
  }

  componentDidMount() {
    const { svg } = setupScatterplot(width, height, margin);
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

    window.addEventListener("resize", this.handleResize);
    const viewport = this.getViewportDimensions();

    this.setState({
      regl,
      flagBuffer,
      pointBuffer,
      colorBuffer,
      svg,
      drawPoints,
      projectionTF,
      viewport,
    });

    this.updateState(null);
  }

  componentDidUpdate(prevProps) {
    this.updateState(prevProps);
  }

  componentWillUnmount() {
    window.removeEventListener("resize", this.updateViewportDimensions);
  }

  getViewportDimensions = () => {
    return {
      viewport: {
        height: window.height,
        width: window.width,
      },
    };
  };

  static setupScales(expressionX, expressionY) {
    const { min: xMin, max: xMax } = expressionX.icol(0).summarize();
    const { min: yMin, max: yMax } = expressionY.icol(0).summarize();
    const xScale = d3.scaleLinear().domain([xMin, xMax]).range([0, width]);
    const yScale = d3.scaleLinear().domain([yMin, yMax]).range([height, 0]);

    return {
      xScale,
      yScale,
    };
  }

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

  updateReglState(scatterplotState, colorState, pointDilationState) {
    const { renderCache, state } = this;
    const { svg, regl, pointBuffer, colorBuffer, flagBuffer } = state;

    if (!regl || !svg) return;

    const {
      scatterplotYYaccessor,
      scatterplotXXaccessor,
      expressionXDf,
      expressionYDf,
    } = scatterplotState;
    const xCol = expressionXDf.icol(0);
    const yCol = expressionYDf.icol(0);
    let { xScale, yScale } = renderCache;
    if (
      xCol !== state.scatterplotState.expressionXDf?.icol(0) ||
      yCol !== state.scatterplotState.expressionYDf?.icol(0) ||
      !xScale ||
      !yScale
    ) {
      ({ xScale, yScale } = Scatterplot.setupScales(
        expressionXDf,
        expressionYDf
      ));
      Scatterplot.drawAxesSVG(
        xScale,
        yScale,
        svg,
        scatterplotYYaccessor,
        scatterplotXXaccessor
      );
      renderCache.xScale = xScale;
      renderCache.yScale = yScale;
    }

    const positions = this.computePointPositions(
      xCol.asArray(),
      yCol.asArray(),
      xScale,
      yScale
    );
    if (positions !== renderCache.positions) {
      renderCache.positions = positions;
      pointBuffer({ data: positions, dimension: 2 });
    }

    /* colors for each point */
    const { colors, colorTable, colorDf } = colorState;
    const { colorAccessor } = colors;
    const _colors = this.computePointColors(colorTable.rgb);
    if (_colors !== renderCache.colors) {
      renderCache.colors = _colors;
      colorBuffer({ data: _colors, dimension: 3 });
    }

    /* flags */
    const { pointDilation, pointDilationDf } = pointDilationState;
    const { crossfilter } = this.props;
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
    if (flags !== renderCache.flags) {
      renderCache.flags = flags;
      flagBuffer({ data: flags, dimension: 1 });
    }
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

  async updateState(prevProps) {
    const {
      annoMatrix,
      scatterplotXXaccessor,
      scatterplotYYaccessor,
      colors,
      crossfilter,
      pointDilation,
    } = this.props;
    if (!annoMatrix) return;

    if (
      annoMatrix !== prevProps?.annoMatrix ||
      scatterplotXXaccessor !== prevProps?.scatterplotXXaccessor ||
      scatterplotYYaccessor !== prevProps?.scatterplotYYaccessor ||
      colors !== prevProps?.colors ||
      pointDilation !== prevProps?.pointDilation
    ) {
      this.setState({ status: "pending" });
      try {
        const [
          expressionXDf,
          expressionYDf,
          colorDf,
          pointDilationDf,
        ] = await this.fetchData(
          scatterplotXXaccessor,
          scatterplotYYaccessor,
          colors,
          pointDilation
        );
        const colorTable = this.updateColorTable(colors, colorDf);
        const scatterplotState = {
          scatterplotXXaccessor,
          scatterplotYYaccessor,
          expressionXDf,
          expressionYDf,
        };
        const colorState = { colors, colorDf, colorTable };
        const pointDilationState = { pointDilation, pointDilationDf };
        this.updateReglState(scatterplotState, colorState, pointDilationState);
        this.setState({
          status: "success",
          scatterplotState,
          colorState,
          pointDilationState,
        });
      } catch (error) {
        this.setState({ status: "error" });
        throw error;
      }
      return;
    }

    if (crossfilter !== prevProps?.crossfilter) {
      const { scatterplotState, colorState, pointDilationState } = this.state;
      this.updateReglState(scatterplotState, colorState, pointDilationState);
    }
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

  render() {
    const { dispatch } = this.props;
    const { minimized, status, regl } = this.state;

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
            ref={(canvas) => {
              this.reglCanvas = canvas;
            }}
          />
        </div>
      </div>
    );
  }
}

export default Scatterplot;
