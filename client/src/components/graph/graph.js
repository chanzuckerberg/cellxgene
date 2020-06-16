// jshint esversion: 6
import React from "react";
import * as d3 from "d3";
import { connect } from "react-redux";
import { mat3, vec2 } from "gl-matrix";
import _regl from "regl";
import memoize from "memoize-one";

import setupSVGandBrushElements from "./setupSVGandBrush";
import _camera from "../../util/camera";
import _drawPoints from "./drawPointsRegl";
import { isTypedArray } from "../../util/typeHelpers";
import {
  createColorTable,
  createColorQuery,
} from "../../util/stateManager/colorHelpers";

import GraphOverlayLayer from "./overlays/graphOverlayLayer";
import CentroidLabels from "./overlays/centroidLabels";
import actions from "../../actions";

/*
Simple 2D transforms control all point painting.  There are three:
  * model - convert from underlying per-point coordinate to a layout.
    Currently used to move from data to webgl coordinate system.
  * camera - apply a 2D camera transformation (pan, zoom)
  * projection - apply any transformation required for screen size and layout
*/

function createProjectionTF(viewportWidth, viewportHeight) {
  /*
  the projection transform accounts for the screen size & other layout
  */
  const fractionToUse = 0.95; // fraction of min dimension to use
  const topGutterSizePx = 32; // toolbar box height
  const heightMinusGutter = viewportHeight - topGutterSizePx;
  const minDim = Math.min(viewportWidth, heightMinusGutter);
  const aspectScale = [
    (fractionToUse * minDim) / viewportWidth,
    (fractionToUse * minDim) / viewportHeight,
  ];
  const m = mat3.create();
  mat3.fromTranslation(m, [
    0,
    -topGutterSizePx / viewportHeight / aspectScale[1],
  ]);
  mat3.scale(m, m, aspectScale);
  return m;
}

function createModelTF() {
  /*
  preallocate coordinate system transformation between data and gl.
  Data arrives in a [0,1] range, and we operate elsewhere in [-1,1].
  */
  const m = mat3.fromScaling(mat3.create(), [2, 2]);
  mat3.translate(m, m, [-0.5, -0.5]);
  return m;
}

function renderThrottle(callback) {
  /*
  This wraps a call to requestAnimationFrame(), enforcing a single
  render callback at any given time (ie, you can call this any number
  of times, and it will coallesce multiple inter-frame calls into a
  single render).
  */
  let rafCurrentlyInProgress = null;
  return function f() {
    if (rafCurrentlyInProgress) return;
    const context = this;
    rafCurrentlyInProgress = window.requestAnimationFrame(() => {
      callback.apply(context);
      rafCurrentlyInProgress = null;
    });
  };
}

@connect((state) => ({
  annoMatrix: state.annoMatrix,
  crossfilter: state.obsCrossfilter,
  selectionTool: state.graphSelection.tool,
  currentSelection: state.graphSelection.selection,
  layoutChoice: state.layoutChoice,
  graphInteractionMode: state.controls.graphInteractionMode,
  colors: state.colors,
  pointDilation: state.pointDilation,
}))
class Graph extends React.Component {
  computePointPositions = memoize((X, Y, modelTF) => {
    /*
    compute the model coordinate for each point
    */
    const positions = new Float32Array(2 * X.length);
    for (let i = 0, len = X.length; i < len; i += 1) {
      const p = vec2.fromValues(X[i], Y[i]);
      vec2.transformMat3(p, p, modelTF);
      positions[2 * i] = p[0];
      positions[2 * i + 1] = p[1];
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
    // (world, crossfilter, colorAccessor, pointDilation) => {
    (crossfilter, colorDf, colorAccessor, pointDilationDf, pointDilation) => {
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

      const flagSelected = 1;
      const flagNaN = 2;
      const flagHighlight = 4;

      const flags = this.computeSelectedFlags(
        crossfilter,
        flagSelected,
        0
      ).slice();

      const colorByData = colorDf?.col(colorAccessor)?.asArray();

      const {
        metadataField: pointDilationCategory,
        categoryField: pointDilationLabel,
      } = pointDilation;
      const highlightData = pointDilationDf
        ?.col(pointDilationCategory)
        ?.asArray();

      if (colorByData || highlightData) {
        for (let i = 0, len = flags.length; i < len; i += 1) {
          if (highlightData) {
            flags[i] +=
              highlightData[i] === pointDilationLabel ? flagHighlight : 0;
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
    const viewport = this.getViewportDimensions();
    this.reglCanvas = null;
    const modelTF = createModelTF();
    this.state = {
      toolSVG: null,
      tool: null,
      container: null,
      viewport,

      // projection
      camera: null,
      modelTF,
      modelInvTF: mat3.invert([], modelTF),

      projectionTF: null,

      // regl state
      regl: null,
      drawPoints: null,
      pointBuffer: null,
      colorBuffer: null,
      flagBuffer: null,

      // dataframes
      layoutDf: null,
      colorDf: null,
      pointDilationDf: null,

      // color lookup table
      colorTable: null,

      // regl buffer data
      positions: null,
      colors: null,
      flags: null,
    };
  }

  updateColorTable(colorDf) {
    /* update color table state */
    const { annoMatrix, colors } = this.props;
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

  updateReglState(layoutDf, colorDf, colorTable, pointDilationDf) {
    /* update all regl-related state */
    const newState = {};
    let needsRepaint = false;
    const { layoutChoice } = this.props;
    const { colorAccessor } = this.props.colors;
    const { state } = this;
    const { regl, pointBuffer, colorBuffer, flagBuffer, modelTF } = state;

    // still initializing?
    if (!regl) return newState;

    /* point coordinates */
    const X = layoutDf.col(layoutChoice.currentDimNames[0]).asArray();
    const Y = layoutDf.col(layoutChoice.currentDimNames[1]).asArray();
    const positions = this.computePointPositions(X, Y, modelTF);
    if (positions !== state.positions) {
      newState.positions = positions;
      pointBuffer({ data: positions, dimension: 2 });
      needsRepaint = true;
    }

    /* point colors */
    const colors = this.computePointColors(colorTable.rgb);
    if (colors !== state.colors) {
      /* update our cache & GL if the buffer changes */
      newState.colors = colors;
      colorBuffer({ data: colors, dimension: 3 });
      needsRepaint = true;
    }

    /* flags */
    const { crossfilter, pointDilation } = this.props;
    const flags = this.computePointFlags(
      crossfilter,
      colorDf,
      colorAccessor,
      pointDilationDf,
      pointDilation
    );
    if (flags !== state.flags) {
      newState.flags = flags;
      flagBuffer({ data: flags, dimension: 1 });
      needsRepaint = true;
    }

    if (needsRepaint) {
      const { drawPoints, camera, projectionTF } = state;
      this.renderPoints(
        regl,
        drawPoints,
        colorBuffer,
        pointBuffer,
        flagBuffer,
        camera,
        projectionTF
      );
    }

    return newState;
  }

  colorByQuery() {
    const { annoMatrix, colors } = this.props;
    const { schema } = annoMatrix;
    const { colorMode, colorAccessor } = colors;
    return createColorQuery(colorMode, colorAccessor, schema);
  }

  async fetchData() {
    /*
    fetch all data needed.  Includes:
      - the color by dataframe
      - the layout dataframe
      - the point dilation dataframe
    */
    const { annoMatrix, layoutChoice, colors, pointDilation } = this.props;
    const { colorAccessor } = colors;
    const { metadataField: pointDilationAccessor } = pointDilation;

    const promises = [];
    // layout
    promises.push(annoMatrix.fetch("emb", layoutChoice.current));

    // color
    const query = this.colorByQuery();
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
      colors,
      layoutChoice,
      crossfilter,
      pointDilation,
    } = this.props;
    if (!annoMatrix) return;

    if (
      annoMatrix !== prevProps?.annoMatrix ||
      layoutChoice !== prevProps?.layoutChoice ||
      colors !== prevProps?.colors ||
      pointDilation !== prevProps?.pointDilation
    ) {
      this.setState({ status: "pending" });
      try {
        const [layoutDf, colorDf, pointDilationDf] = await this.fetchData();
        const colorTable = this.updateColorTable(colorDf);
        this.setState({
          status: "success",
          layoutDf,
          colorDf,
          colorTable,
          pointDilationDf,
          ...this.updateReglState(
            layoutDf,
            colorDf,
            colorTable,
            pointDilationDf
          ),
        });
      } catch (error) {
        this.setState({ status: "error", error });
        throw error;
      }
      return;
    }

    if (
      crossfilter !== prevProps?.crossfilter ||
      pointDilation !== prevProps?.pointDilation
    ) {
      const { layoutDf, colorDf, colorTable, pointDilationDf } = this.state;
      this.setState(
        this.updateReglState(layoutDf, colorDf, colorTable, pointDilationDf)
      );
    }
  }

  componentDidMount() {
    window.addEventListener("resize", this.handleResize);

    // setup canvas, webgl draw function and camera
    const camera = _camera(this.reglCanvas);
    const regl = _regl(this.reglCanvas);
    const drawPoints = _drawPoints(regl);

    // preallocate webgl buffers
    const pointBuffer = regl.buffer();
    const colorBuffer = regl.buffer();
    const flagBuffer = regl.buffer();

    // create all default rendering transformations
    const projectionTF = createProjectionTF(
      this.reglCanvas.width,
      this.reglCanvas.height
    );

    this.setState({
      regl,
      drawPoints,
      pointBuffer,
      colorBuffer,
      flagBuffer,
      camera,
      projectionTF,
    });

    this.updateState(null);
  }

  componentDidUpdate(prevProps, prevState) {
    const {
      crossfilter,
      selectionTool,
      currentSelection,
      layoutChoice,
      graphInteractionMode,
      pointDilation,
    } = this.props;
    const { regl, toolSVG, camera, modelTF, viewport } = this.state;
    let { projectionTF } = this.state;
    const { reglCanvas } = this;
    const hasResized =
      reglCanvas &&
      (prevState.viewport.height !== reglCanvas.height ||
        prevState.viewport.width !== reglCanvas.width);
    let stateChanges = {};
    let needsRepaint = hasResized;

    if (hasResized) {
      projectionTF = createProjectionTF(reglCanvas.width, reglCanvas.height);
      stateChanges = {
        ...stateChanges,
        projectionTF,
      };
    }

    this.updateState(prevProps);

    if (hasResized) {
      // If the window size has changed we want to recreate all SVGs
      stateChanges = {
        ...stateChanges,
        ...this.createToolSVG(),
      };
    } else if (
      (viewport.height && viewport.width && !toolSVG) ||
      selectionTool !== prevProps.selectionTool
    ) {
      // first time or change of selection tool
      stateChanges = { ...stateChanges, ...this.createToolSVG() };
    } else if (prevProps.graphInteractionMode !== graphInteractionMode) {
      // If lasso/zoom is switched
      stateChanges = {
        ...stateChanges,
        ...this.createToolSVG(),
      };
    }

    /*
    if the selection tool or state has changed, ensure that the selection
    tool correctly reflects the underlying selection.
    */
    if (
      currentSelection !== prevProps.currentSelection ||
      graphInteractionMode !== prevProps.graphInteractionMode ||
      stateChanges.toolSVG
    ) {
      const { tool, container } = this.state;
      this.selectionToolUpdate(
        stateChanges.tool ? stateChanges.tool : tool,
        stateChanges.container ? stateChanges.container : container
      );
    }
    if (Object.keys(stateChanges).length > 0) {
      // Preventing update loop via stateChanges and diff checks
      // eslint-disable-next-line react/no-did-update-set-state
      this.setState(stateChanges);
    }
  }

  componentWillUnmount() {
    window.removeEventListener("resize", this.handleResize);
  }

  handleResize = () => {
    const { state } = this.state;
    const viewport = this.getViewportDimensions();
    this.setState({
      ...state,
      viewport,
    });
  };

  getViewportDimensions = () => {
    const { viewportRef } = this.props;
    return {
      height: viewportRef.clientHeight,
      width: viewportRef.clientWidth,
    };
  };

  handleCanvasEvent = (e) => {
    const { camera, projectionTF } = this.state;
    if (e.type !== "wheel") e.preventDefault();
    if (camera.handleEvent(e, projectionTF)) {
      this.renderCanvas();
      this.setState((state) => {
        return { ...state, updateOverlay: !state.updateOverlay };
      });
    }
  };

  createToolSVG = () => {
    /*
    Called from componentDidUpdate. Create the tool SVG, and return any
    state changes that should be passed to setState().
    */
    const { selectionTool, graphInteractionMode } = this.props;
    const { viewport } = this.state;

    /* clear out whatever was on the div, even if nothing, but usually the brushes etc */

    d3.select("#lasso-layer").selectAll(".lasso-group").remove();

    // Don't render or recreate toolSVG if currently in zoom mode
    if (graphInteractionMode !== "select") {
      // don't return "change" of state unless we are really changing it!
      const { toolSVG } = this.state;
      if (toolSVG === undefined) return {};
      return { toolSVG: undefined };
    }

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

    const { svg: newToolSVG, tool, container } = setupSVGandBrushElements(
      selectionTool,
      handleStart,
      handleDrag,
      handleEnd,
      handleCancel,
      viewport
    );

    return { toolSVG: newToolSVG, tool, container };
  };

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
          this.mapPointToScreen(currentSelection.brushCoords.southeast),
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

  lassoToolUpdate(tool) {
    /*
    this is called from componentDidUpdate(), so be very careful using
    anything from this.state, which may be updated asynchronously.
    */
    const { currentSelection } = this.props;
    if (currentSelection.mode === "within-polygon") {
      /*
      if there is a current selection, make sure the lasso tool matches
      */
      const polygon = currentSelection.polygon.map((p) =>
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

  mapScreenToPoint(pin) {
    /*
    Map an XY coordinates from screen domain to cell/point range,
    accounting for current pan/zoom camera.
    */

    const { camera, projectionTF, modelInvTF, viewport } = this.state;
    const cameraInvTF = camera.invView();

    /* screen -> gl */
    const x = (2 * pin[0]) / viewport.width - 1;
    const y = 2 * (1 - pin[1] / viewport.height) - 1;

    const xy = vec2.fromValues(x, y);
    const projectionInvTF = mat3.invert(mat3.create(), projectionTF);
    vec2.transformMat3(xy, xy, projectionInvTF);
    vec2.transformMat3(xy, xy, cameraInvTF);
    vec2.transformMat3(xy, xy, modelInvTF);
    return xy;
  }

  mapPointToScreen(xyCell) {
    /*
    Map an XY coordinate from cell/point domain to screen range.  Inverse
    of mapScreenToPoint()
    */

    const { camera, projectionTF, modelTF, viewport } = this.state;
    const cameraTF = camera.view();

    const xy = vec2.transformMat3(vec2.create(), xyCell, modelTF);
    vec2.transformMat3(xy, xy, cameraTF);
    vec2.transformMat3(xy, xy, projectionTF);

    return [
      Math.round(((xy[0] + 1) * viewport.width) / 2),
      Math.round(-((xy[1] + 1) / 2 - 1) * viewport.height),
    ];
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

    const { dispatch, layoutChoice } = this.props;
    const s = d3.event.selection;
    const northwest = this.mapScreenToPoint(s[0]);
    const southeast = this.mapScreenToPoint(s[1]);
    const [minX, maxY] = northwest;
    const [maxX, minY] = southeast;
    dispatch(
      actions.graphBrushChangeAction(layoutChoice.current, {
        minX,
        minY,
        maxX,
        maxY,
        northwest,
        southeast,
      })
    );
  }

  handleBrushStartAction() {
    // Ignore programatically generated events.
    if (!d3.event.sourceEvent) return;

    const { dispatch } = this.props;
    dispatch(actions.graphBrushStartAction());
  }

  handleBrushEndAction() {
    // Ignore programatically generated events.
    if (!d3.event.sourceEvent) return;

    /*
    coordinates will be included if selection made, null
    if selection cleared.
    */
    const { dispatch, layoutChoice } = this.props;
    const s = d3.event.selection;
    if (s) {
      const northwest = this.mapScreenToPoint(s[0]);
      const southeast = this.mapScreenToPoint(s[1]);
      const [minX, maxY] = northwest;
      const [maxX, minY] = southeast;
      dispatch(
        actions.graphBrushEndAction(layoutChoice.current, {
          minX,
          minY,
          maxX,
          maxY,
          northwest,
          southeast,
        })
      );
    } else {
      dispatch(actions.graphBrushDeselectAction(layoutChoice.current));
    }
  }

  handleBrushDeselectAction() {
    const { dispatch, layoutChoice } = this.props;
    dispatch(actions.graphBrushDeselectAction(layoutChoice.current));
  }

  handleLassoStart() {
    const { dispatch, layoutChoice } = this.props;
    dispatch(actions.graphLassoStartAction(layoutChoice.current));
  }

  // when a lasso is completed, filter to the points within the lasso polygon
  handleLassoEnd(polygon) {
    const minimumPolygonArea = 10;
    const { dispatch, layoutChoice } = this.props;

    if (
      polygon.length < 3 ||
      Math.abs(d3.polygonArea(polygon)) < minimumPolygonArea
    ) {
      // if less than three points, or super small area, treat as a clear selection.
      dispatch(actions.graphLassoDeselectAction(layoutChoice.current));
    } else {
      dispatch(
        actions.graphLassoEndAction(
          layoutChoice.current,
          polygon.map((xy) => this.mapScreenToPoint(xy))
        )
      );
    }
  }

  handleLassoCancel() {
    const { dispatch, layoutChoice } = this.props;
    dispatch(actions.graphLassoCancelAction(layoutChoice.current));
  }

  handleLassoDeselectAction() {
    const { dispatch, layoutChoice } = this.props;
    dispatch(actions.graphLassoDeselectAction(layoutChoice.current));
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
      data: e.target.value,
    });
  }

  renderPoints(
    regl,
    drawPoints,
    colorBuffer,
    pointBuffer,
    flagBuffer,
    camera,
    projectionTF
  ) {
    const { annoMatrix } = this.props;
    if (!this.reglCanvas || !annoMatrix) return;

    const { schema } = annoMatrix;
    const cameraTF = camera.view();
    const projView = mat3.multiply(mat3.create(), projectionTF, cameraTF);
    const { width, height } = this.reglCanvas;
    regl.poll();
    regl.clear({
      depth: 1,
      color: [1, 1, 1, 1],
    });
    drawPoints({
      distance: camera.distance(),
      color: colorBuffer,
      position: pointBuffer,
      flag: flagBuffer,
      count: annoMatrix.nObs,
      projView,
      nPoints: schema.dataframe.nObs,
      minViewportDimension: Math.min(width, height),
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
      camera,
      projectionTF,
    } = this.state;
    this.renderPoints(
      regl,
      drawPoints,
      colorBuffer,
      pointBuffer,
      flagBuffer,
      camera,
      projectionTF
    );
  });

  render() {
    const { graphInteractionMode } = this.props;
    const { modelTF, projectionTF, camera, viewport } = this.state;
    const cameraTF = camera?.view()?.slice();

    const { status } = this.state;
    if (status === "error") return null;

    return (
      <div
        id="graph-wrapper"
        style={{
          position: "relative",
          top: 0,
          left: 0,
        }}
      >
        <GraphOverlayLayer
          width={viewport.width}
          height={viewport.height}
          cameraTF={cameraTF}
          modelTF={modelTF}
          projectionTF={projectionTF}
          handleCanvasEvent={
            graphInteractionMode === "zoom" ? this.handleCanvasEvent : undefined
          }
        >
          <CentroidLabels />
        </GraphOverlayLayer>
        <svg
          id="lasso-layer"
          data-testid="layout-overlay"
          className="graph-svg"
          style={{
            position: "absolute",
            top: 0,
            left: 0,
            zIndex: 1,
          }}
          width={viewport.width}
          height={viewport.height}
          pointerEvents={graphInteractionMode === "select" ? "auto" : "none"}
        />
        <canvas
          width={viewport.width}
          height={viewport.height}
          style={{
            position: "absolute",
            top: 0,
            left: 0,
            padding: 0,
            margin: 0,
            shapeRendering: "crispEdges",
          }}
          className="graph-canvas"
          data-testid="layout-graph"
          ref={(canvas) => {
            this.reglCanvas = canvas;
          }}
          onMouseDown={this.handleCanvasEvent}
          onMouseUp={this.handleCanvasEvent}
          onMouseMove={this.handleCanvasEvent}
          onDoubleClick={this.handleCanvasEvent}
          onWheel={this.handleCanvasEvent}
        />
      </div>
    );
  }
}

export default Graph;
