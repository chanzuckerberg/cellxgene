// jshint esversion: 6
import React from "react";
import * as d3 from "d3";
import { connect } from "react-redux";
import mat4 from "gl-mat4";
import _regl from "regl";
import memoize from "memoize-one";

import * as globals from "../../globals";
import setupSVGandBrushElements from "./setupSVGandBrush";
import _camera from "../../util/camera";
import _drawPoints from "./drawPointsRegl";
import scaleLinear from "../../util/scaleLinear";

/* https://bl.ocks.org/mbostock/9078690 - quadtree for onClick / hover selections */

@connect(state => ({
  world: state.world,
  universe: state.universe,
  crossfilter: state.crossfilter,
  responsive: state.responsive,
  colorRGB: state.colors.rgb,
  opacityForDeselectedCells: state.controls.opacityForDeselectedCells,
  userDefinedGenes: state.controls.userDefinedGenes,
  diffexpGenes: state.controls.diffexpGenes,
  colorAccessor: state.colors.colorAccessor,
  scatterplotXXaccessor: state.controls.scatterplotXXaccessor,
  scatterplotYYaccessor: state.controls.scatterplotYYaccessor,
  celllist1: state.differential.celllist1,
  celllist2: state.differential.celllist2,
  libraryVersions: state.config?.library_versions, // eslint-disable-line camelcase
  undoDisabled: state["@@undoable/past"].length === 0,
  redoDisabled: state["@@undoable/future"].length === 0,
  selectionTool: state.graphSelection.tool,
  currentSelection: state.graphSelection.selection,
  layoutChoice: state.layoutChoice,
  graphInteractionMode: state.controls.graphInteractionMode
}))
class Graph extends React.Component {
  computePointPositions = memoize((X, Y, scaleX, scaleY) => {
    /*
    compute webgl coordinate buffer for each point
    */
    const positions = new Float32Array(2 * X.length);
    for (let i = 0, len = X.length; i < len; i += 1) {
      positions[2 * i] = scaleX(X[i]);
      positions[2 * i + 1] = scaleY(Y[i]);
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

  computePointSizes = memoize((len, crossfilter) => {
    /*
    compute webgl dot size for each point
    */
    const sizes = new Float32Array(len);
    crossfilter.fillByIsSelected(sizes, 4, 0.2);
    return sizes;
  });

  constructor(props) {
    super(props);
    this.count = 0;
    this.graphPaddingTop = 80;
    this.graphPaddingBottom = 45;
    this.graphPaddingRight = globals.leftSidebarWidth;
    this.renderCache = {
      X: null,
      Y: null,
      positions: null,
      colors: null,
      sizes: null
    };
    this.state = {
      svg: null,
      tool: null,
      container: null
    };
  }

  componentDidMount() {
    // setup canvas and camera
    const camera = _camera(this.reglCanvas, { scale: true, rotate: false });
    const regl = _regl(this.reglCanvas);

    const drawPoints = _drawPoints(regl);

    // preallocate buffers
    const pointBuffer = regl.buffer();
    const colorBuffer = regl.buffer();
    const sizeBuffer = regl.buffer();

    // preallocate coordinate system transformation between data and gl
    const fractionToUse = 0.98; // fraction of dimension to use
    const transform = {
      glScaleX: scaleLinear([0, 1], [-1 * fractionToUse, 1 * fractionToUse]),
      glScaleY: scaleLinear([0, 1], [1 * fractionToUse, -1 * fractionToUse])
    };

    /* first time, but this duplicates above function, should be possile to avoid this */
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
      drawPoints,
      pointBuffer,
      colorBuffer,
      sizeBuffer,
      camera,
      reglRender,
      transform
    });
  }

  componentDidUpdate(prevProps) {
    const { renderCache } = this;
    const {
      world,
      crossfilter,
      colorRGB,
      responsive,
      selectionTool,
      currentSelection,
      layoutChoice,
      graphInteractionMode
    } = this.props;
    const { reglRender, regl, svg } = this.state;
    let stateChanges = {};

    if (
      reglRender &&
      this.reglRenderState === "rendering" &&
      graphInteractionMode !== "zoom"
    ) {
      reglRender.cancel();
      this.reglRenderState = "paused";
    }

    if (
      reglRender &&
      this.reglRenderState !== "rendering" &&
      graphInteractionMode === "zoom"
    ) {
      this.restartReglLoop();
      this.reglRenderState = "rendering";
    }

    if (regl && world) {
      /* update the regl state */
      const { obsLayout, nObs } = world;
      const {
        drawPoints,
        transform,
        camera,
        pointBuffer,
        colorBuffer,
        sizeBuffer
      } = this.state;

      /* coordinates for each point */
      const { glScaleX, glScaleY } = transform;
      const X = obsLayout.col(layoutChoice.currentDimNames[0]).asArray();
      const Y = obsLayout.col(layoutChoice.currentDimNames[1]).asArray();
      const newPositions = this.computePointPositions(X, Y, glScaleX, glScaleY);
      if (renderCache.positions !== newPositions) {
        /* update our cache & GL if the buffer changes */
        renderCache.positions = newPositions;
        pointBuffer({ data: newPositions, dimension: 2 });
      }

      /* colors for each point */
      const newColors = this.computePointColors(colorRGB);
      if (renderCache.colors !== newColors) {
        /* update our cache & GL if the buffer changes */
        renderCache.colors = newColors;
        colorBuffer({ data: newColors, dimension: 3 });
      }

      /* sizes for each point */
      const newSizes = this.computePointSizes(nObs, crossfilter);
      if (renderCache.sizes !== newSizes) {
        /* update our cache & GL if the buffer changes */
        renderCache.size = newSizes;
        sizeBuffer({ data: newSizes, dimension: 1 });
      }

      this.count = nObs;

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

    if (
      prevProps.responsive.height !== responsive.height ||
      prevProps.responsive.width !== responsive.width ||
      /* first time */
      (responsive.height && responsive.width && !svg) ||
      selectionTool !== prevProps.selectionTool
    ) {
      /* clear out whatever was on the div, even if nothing, but usually the brushes etc */
      d3.select("#graphAttachPoint")
        .selectAll("svg")
        .remove();

      let handleStart;
      let handleDrag;
      let handleEnd;
      let handleCancel;
      if (selectionTool === "brush") {
        handleStart = this.handleBrushStartAction.bind(this);
        handleDrag = this.handleBrushDragAction.bind(this);
        handleEnd = this.handleBrushEndAction.bind(this);
      } else {
        handleStart = this.handleLassoStart.bind(this);
        handleEnd = this.handleLassoEnd.bind(this);
        handleCancel = this.handleLassoCancel.bind(this);
      }
      const { svg: newSvg, tool, container } = setupSVGandBrushElements(
        selectionTool,
        handleStart,
        handleDrag,
        handleEnd,
        handleCancel,
        responsive,
        this.graphPaddingRight
      );
      stateChanges = { ...stateChanges, svg: newSvg, tool, container };
    }

    /*
    if the selection tool or state has changed, ensure that the selection
    tool correctly reflects the underlying selection.
    */
    if (
      currentSelection !== prevProps.currentSelection ||
      graphInteractionMode !== prevProps.graphInteractionMode ||
      stateChanges.svg
    ) {
      const { tool, container } = this.state;
      this.selectionToolUpdate(
        stateChanges.tool ? stateChanges.tool : tool,
        stateChanges.container ? stateChanges.container : container
      );
    }

    if (Object.keys(stateChanges).length > 0) {
      this.setState(stateChanges);
    }
  }

  brushToolUpdate(tool, container) {
    /*
    this is called from componentDidUpdate(), so be very careful using
    anything from this.state, which may be updated asynchronously.
    */
    const { currentSelection } = this.props;
    if (container) {
      const toolCurrentSelection = d3.brushSelection(container.node());

      if (currentSelection.mode === "within-rect") {
        /*
        if there is a selection, make sure the brush tool matches
        */
        const screenCoords = [
          this.mapPointToScreen(currentSelection.brushCoords.northwest),
          this.mapPointToScreen(currentSelection.brushCoords.southeast)
        ];
        if (!toolCurrentSelection) {
          /* tool is not selected, so just move the brush */
          container.call(tool.move, screenCoords);
        } else {
          /* there is an active selection and a brush - make sure they match */
          /* this just sums the difference of each dimension, of each point */
          let delta = 0;
          for (let x = 0; x < 2; x += 1) {
            for (let y = 0; y < 2; y += 1) {
              delta += Math.abs(
                screenCoords[x][y] - toolCurrentSelection[x][y]
              );
            }
          }
          if (delta > 0) {
            container.call(tool.move, screenCoords);
          }
        }
      } else if (toolCurrentSelection) {
        /* no selection, so clear the brush tool if it is set */
        container.call(tool.move, null);
      }
    }
  }

  lassoToolUpdate(tool, container) {
    /*
    this is called from componentDidUpdate(), so be very careful using
    anything from this.state, which may be updated asynchronously.
    */
    const { currentSelection } = this.props;
    if (currentSelection.mode === "within-polygon") {
      /*
      if there is a current selection, make sure the lasso tool matches
      */
      const polygon = currentSelection.polygon.map(p =>
        this.mapPointToScreen(p)
      );
      tool.move(polygon);
    } else {
      tool.reset();
    }
  }

  selectionToolUpdate(tool, container) {
    /*
    this is called from componentDidUpdate(), so be very careful using
    anything from this.state, which may be updated asynchronously.
    */
    const { selectionTool } = this.props;
    switch (selectionTool) {
      case "brush":
        this.brushToolUpdate(tool, container);
        break;
      case "lasso":
        this.lassoToolUpdate(tool, container);
        break;
      default:
        /* punt? */
        break;
    }
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

  restartReglLoop() {
    const {
      regl,
      drawPoints,
      sizeBuffer,
      colorBuffer,
      pointBuffer,
      camera
    } = this.state;
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
      reglRender
    });
  }

  mapScreenToPoint(pin) {
    /*
    Map an XY coordinates from screen domain to cell/point range,
    accounting for current pan/zoom camera.
    */
    const { responsive } = this.props;
    const { regl, camera, transform } = this.state;
    const { glScaleX, glScaleY } = transform;

    const gl = regl._gl;

    // get aspect ratio
    const aspect = gl.drawingBufferWidth / gl.drawingBufferHeight;
    const scale = aspect < 1 ? 1 / aspect : 1;

    // compute inverse view matrix
    const inverse = mat4.invert([], camera.view());

    // transform screen coordinates -> cell coordinates
    const x = (2 * pin[0]) / (responsive.width - this.graphPaddingRight) - 1;
    const y = 2 * (1 - pin[1] / (responsive.height - this.graphPaddingTop)) - 1;
    const pout = [
      x * inverse[14] * aspect * scale + inverse[12],
      -(y * inverse[14] * scale + inverse[13])
    ];

    const xy = [glScaleX.invert(pout[0]), glScaleY.invert(pout[1])];
    return xy;
  }

  mapPointToScreen(xyCell) {
    /*
    Map an XY coordinate from cell/point domain to screen range.  Inverse
    of mapScreenToPoint()
    */
    const { responsive } = this.props;
    const { regl, camera, transform } = this.state;
    const { glScaleX, glScaleY } = transform;

    const gl = regl._gl;

    // get aspect ratio
    const aspect = gl.drawingBufferWidth / gl.drawingBufferHeight;
    const scale = aspect < 1 ? 1 / aspect : 1;

    // compute inverse view matrix
    let inverse = mat4.invert([], camera.view());

    // variable names are choosen to reflect inverse of those used
    // in mapScreenToPoint().
    const pout = [glScaleX(xyCell[0]), glScaleY(xyCell[1])];
    const x = (pout[0] - inverse[12]) / aspect / scale / inverse[14];
    const y = (-pout[1] - inverse[13]) / scale / inverse[14];

    const pin = [
      Math.round(((x + 1) * (responsive.width - this.graphPaddingRight)) / 2),
      Math.round(
        -((y + 1) / 2 - 1) * (responsive.height - this.graphPaddingTop)
      )
    ];
    return pin;
  }

  handleBrushDragAction() {
    /*
      event describing brush position:
      @-------|
      |       |
      |       |
      |-------@
    */
    // ignore programatically generated events
    if (d3.event.sourceEvent === null || !d3.event.selection) return;

    const { dispatch } = this.props;
    const s = d3.event.selection;
    const brushCoords = {
      northwest: this.mapScreenToPoint([s[0][0], s[0][1]]),
      southeast: this.mapScreenToPoint([s[1][0], s[1][1]])
    };

    dispatch({
      type: "graph brush change",
      brushCoords
    });
  }

  handleBrushStartAction() {
    // Ignore programatically generated events.
    if (!d3.event.sourceEvent) return;

    const { dispatch } = this.props;
    dispatch({ type: "graph brush start" });
  }

  handleBrushEndAction() {
    // Ignore programatically generated events.
    if (!d3.event.sourceEvent) return;

    /*
    coordinates will be included if selection made, null
    if selection cleared.
    */
    const { dispatch } = this.props;
    const s = d3.event.selection;
    if (s) {
      const brushCoords = {
        northwest: this.mapScreenToPoint(s[0]),
        southeast: this.mapScreenToPoint(s[1])
      };
      dispatch({
        type: "graph brush end",
        brushCoords
      });
    } else {
      dispatch({
        type: "graph brush deselect"
      });
    }
  }

  handleBrushDeselectAction() {
    const { dispatch } = this.props;
    dispatch({
      type: "graph brush deselect"
    });
  }

  handleLassoStart() {
    const { dispatch } = this.props;
    dispatch({
      type: "graph lasso start"
    });
  }

  // when a lasso is completed, filter to the points within the lasso polygon
  handleLassoEnd(polygon) {
    const minimumPolygonArea = 10;
    const { dispatch } = this.props;

    if (
      polygon.length < 3 ||
      Math.abs(d3.polygonArea(polygon)) < minimumPolygonArea
    ) {
      // if less than three points, or super small area, treat as a clear selection.
      dispatch({ type: "graph lasso deselect" });
    } else {
      dispatch({
        type: "graph lasso end",
        polygon: polygon.map(xy => this.mapScreenToPoint(xy)) // transform the polygon
      });
    }
  }

  handleLassoCancel() {
    const { dispatch } = this.props;
    dispatch({ type: "graph lasso cancel" });
  }

  handleLassoDeselectAction() {
    const { dispatch } = this.props;
    dispatch({ type: "graph lasso deselect" });
  }

  handleDeselectAction() {
    const { selectionTool } = this.props;
    if (selectionTool === "brush") this.handleBrushDeselectAction();
    if (selectionTool === "lasso") this.handleLassoDeselectAction();
  }

  handleOpacityRangeChange(e) {
    const { dispatch } = this.props;
    dispatch({
      type: "change opacity deselected cells in 2d graph background",
      data: e.target.value
    });
  }

  render() {
    const { responsive, graphInteractionMode } = this.props;

    return (
      <div id="graphWrapper">
        <div
          style={{
            zIndex: -9999,
            position: "fixed",
            bottom: 0,
            right: 0
          }}
        >
          <div
            style={{
              display: graphInteractionMode === "select" ? "inherit" : "none"
            }}
            id="graphAttachPoint"
          />
          <div style={{ padding: 0, margin: 0 }}>
            <canvas
              width={responsive.width - this.graphPaddingRight}
              height={responsive.height - this.graphPaddingTop}
              data-testid="layout-graph"
              ref={canvas => {
                this.reglCanvas = canvas;
              }}
            />
          </div>
        </div>
      </div>
    );
  }
}

export default Graph;
