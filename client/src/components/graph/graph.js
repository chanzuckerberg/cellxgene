// jshint esversion: 6
import React from "react";
import * as d3 from "d3";
import { connect } from "react-redux";
import { mat3, vec2 } from "gl-matrix";
import _regl from "regl";
import memoize from "memoize-one";

import * as globals from "../../globals";
import setupSVGandBrushElements from "./setupSVGandBrush";
import setupCentroidSVG from "./setupCentroidSVG";
import _camera from "../../util/camera";
import _drawPoints from "./drawPointsRegl";
import { isTypedArray } from "../../util/typeHelpers";

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
    (fractionToUse * minDim) / viewportHeight
  ];
  const m = mat3.create();
  mat3.fromTranslation(m, [
    0,
    -topGutterSizePx / viewportHeight / aspectScale[1]
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

@connect(state => ({
  universe: state.universe,
  world: state.world,
  crossfilter: state.crossfilter,
  responsive: state.responsive,
  colorRGB: state.colors.rgb,
  selectionTool: state.graphSelection.tool,
  currentSelection: state.graphSelection.selection,
  layoutChoice: state.layoutChoice,
  centroidLabel: state.centroidLabel,
  graphInteractionMode: state.controls.graphInteractionMode,
  colorAccessor: state.colors.colorAccessor
}))
class Graph extends React.PureComponent {
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
    this.graphPaddingTop = 0;
    this.graphPaddingRightLeft = globals.leftSidebarWidth * 2;
    this.renderCache = {
      X: null,
      Y: null,
      positions: null,
      colors: null,
      sizes: null,
      flags: null
    };
    this.state = {
      toolSVG: null,
      centroidSVG: null,
      tool: null,
      container: null
    };
  }

  componentDidMount() {
    // setup canvas, webgl draw function and camera
    const camera = _camera(this.reglCanvas);
    const regl = _regl(this.reglCanvas);
    const drawPoints = _drawPoints(regl);

    // preallocate webgl buffers
    const pointBuffer = regl.buffer();
    const colorBuffer = regl.buffer();
    const flagBuffer = regl.buffer();

    // create all default rendering transformations
    const modelTF = createModelTF();
    const projectionTF = createProjectionTF(
      this.reglCanvas.width,
      this.reglCanvas.height
    );

    // initial draw to canvas
    this.renderPoints(
      regl,
      drawPoints,
      colorBuffer,
      pointBuffer,
      flagBuffer,
      camera,
      projectionTF
    );

    this.setState({
      regl,
      drawPoints,
      pointBuffer,
      colorBuffer,
      flagBuffer,
      camera,
      modelTF,
      modelInvTF: mat3.invert([], modelTF),
      projectionTF
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
      graphInteractionMode,
      colorAccessor,
      centroidLabel
    } = this.props;
    const { regl, toolSVG, centroidSVG } = this.state;
    let stateChanges = {};

    if (regl && world) {
      /* update the regl and point rendering state */
      const { obsLayout, nObs } = world;
      const {
        drawPoints,
        camera,
        pointBuffer,
        colorBuffer,
        flagBuffer,
        modelTF
      } = this.state;
      let { projectionTF } = this.state;
      let needsRepaint = false;

      if (
        prevProps.responsive.height !== responsive.height ||
        prevProps.responsive.width !== responsive.width
      ) {
        projectionTF = createProjectionTF(
          this.reglCanvas.width,
          this.reglCanvas.height
        );
        needsRepaint = true;
        stateChanges = {
          ...stateChanges,
          projectionTF
        };
      }

      /* coordinates for each point */
      const X = obsLayout.col(layoutChoice.currentDimNames[0]).asArray();
      const Y = obsLayout.col(layoutChoice.currentDimNames[1]).asArray();
      const newPositions = this.computePointPositions(X, Y, modelTF);
      if (renderCache.positions !== newPositions) {
        /* update our cache & GL if the buffer changes */
        renderCache.positions = newPositions;
        pointBuffer({ data: newPositions, dimension: 2 });
        needsRepaint = true;
      }

      /* colors for each point */
      const newColors = this.computePointColors(colorRGB);
      if (renderCache.colors !== newColors) {
        /* update our cache & GL if the buffer changes */
        renderCache.colors = newColors;
        colorBuffer({ data: newColors, dimension: 3 });
        needsRepaint = true;
      }

      /* flags for each point */
      const newFlags = this.computePointFlags(
        world,
        crossfilter,
        colorAccessor,
        centroidLabel
      );
      if (renderCache.flags !== newFlags) {
        renderCache.flags = newFlags;
        flagBuffer({ data: newFlags, dimension: 1 });
        needsRepaint = true;
      }

      this.count = nObs;

      if (needsRepaint) {
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
    }

    const createToolSVG = () => {
      /* clear out whatever was on the div, even if nothing, but usually the brushes etc */
      d3.select("#graphAttachPoint")
        .select("#tool")
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

      const { svg: newToolSVG, tool, container } = setupSVGandBrushElements(
        selectionTool,
        handleStart,
        handleDrag,
        handleEnd,
        handleCancel,
        responsive,
        this.graphPaddingRight,
        graphInteractionMode
      );

      stateChanges = { ...stateChanges, toolSVG: newToolSVG, tool, container };
    };

    const createCentroidSVG = (
      viewportChange = true /* this should be defaulted to false */
    ) => {
      d3.select("#graphAttachPoint")
        .select("#centroid-container")
        .remove();

      if (
        centroidLabel.labeledCategory === "" ||
        !centroidLabel.labels ||
        graphInteractionMode === "zoom"
      ) {
        return;
      }

      const iter = centroidLabel.labels.entries();
      let pair = iter.next().value;
      let value;
      while (pair) {
        value = pair[1];
        if (!value[2] || viewportChange) {
          console.log("before:", value);

          value[2] = true;

          value.splice(3, 2, ...this.mapPointToScreen([value[0], value[1]]));

          console.log(value);
        }
        pair = iter.next().value;
      }

      const newCentroidSVG = setupCentroidSVG(
        responsive,
        this.graphPaddingRight,
        centroidLabel.labels,
        colorAccessor
      );

      stateChanges = { ...stateChanges, centroidSVG: newCentroidSVG };
    };

    if (
      prevProps.responsive.height !== responsive.height ||
      prevProps.responsive.width !== responsive.width
    ) {
      // If the window size has changed we want to recreate all SVGs
      createToolSVG();
      createCentroidSVG(true);
      // stateChanges = { ...stateChanges, ...this.createToolSVG() };
    } else if (
      (responsive.height && responsive.width && !toolSVG) ||
      selectionTool !== prevProps.selectionTool
    ) {
      // first time or change of selection tool6
      createToolSVG();
    } else if (prevProps.graphInteractionMode !== graphInteractionMode) {
      createToolSVG();
      createCentroidSVG();
    } else if (
      centroidLabel !== prevProps.centroidLabel ||
      (responsive.height && responsive.width && !centroidSVG)
    ) {
      // First time for centroid or label change
      createCentroidSVG();
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
      this.setState(stateChanges);
    }
  }

  handleCanvasEvent = e => {
    const { camera, projectionTF } = this.state;
    if (e.type !== "wheel") e.preventDefault();
    if (camera.handleEvent(e, projectionTF)) {
      this.renderCanvas();
    }
  };

  createToolSVG() {
    /*
    Called from componentDidUpdate. Create the tool SVG, and return any
    state changes that should be passed to setState().
    */
    const { responsive, selectionTool, graphInteractionMode } = this.props;

    /* clear out whatever was on the div, even if nothing, but usually the brushes etc */
    d3.select("#graphAttachPoint")
      .select("#tool")
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

    const { svg: newToolSVG, tool, container } = setupSVGandBrushElements(
      selectionTool,
      handleStart,
      handleDrag,
      handleEnd,
      handleCancel,
      responsive,
      this.graphPaddingRightLeft,
      graphInteractionMode
    );

    return { toolSVG: newToolSVG, tool, container };
  }

  createCentroidSVG() {
    /*
    Called from componentDidUpdate. Create the centroid SVG, and return any
    state changes that should be passed to setState().

    CURRENTLY UNUSED
    */
    const { responsive, centroidLabel, colorAccessor } = this.props;
    d3.select("#graphAttachPoint")
      .select("#centroid-container")
      .remove();

    if (centroidLabel.metadataField === "" || !centroidLabel.centroidXY) {
      return {};
    }

    const centroidScreen = this.mapPointToScreen(centroidLabel.centroidXY);

    const newCentroidSVG = setupCentroidSVG(
      responsive,
      this.graphPaddingRightLeft,
      centroidScreen,
      centroidLabel.categoryField,
      colorAccessor
    );

    return { centroidSVG: newCentroidSVG };
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

  mapScreenToPoint(pin) {
    /*
    Map an XY coordinates from screen domain to cell/point range,
    accounting for current pan/zoom camera.
    */

    const { responsive } = this.props;
    const { camera, projectionTF, modelInvTF } = this.state;
    const cameraInvTF = camera.invView();

    /* screen -> gl */
    const x =
      (2 * pin[0]) / (responsive.width - this.graphPaddingRightLeft) - 1;
    const y = 2 * (1 - pin[1] / (responsive.height - this.graphPaddingTop)) - 1;

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

    const { responsive } = this.props;
    const { camera, projectionTF, modelTF } = this.state;
    const cameraTF = camera.view();

    const xy = vec2.transformMat3(vec2.create(), xyCell, modelTF);
    vec2.transformMat3(xy, xy, cameraTF);
    vec2.transformMat3(xy, xy, projectionTF);

    const pin = [
      Math.round(
        ((xy[0] + 1) * (responsive.width - this.graphPaddingRightLeft)) / 2
      ),
      Math.round(
        -((xy[1] + 1) / 2 - 1) * (responsive.height - this.graphPaddingTop)
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

  renderPoints(
    regl,
    drawPoints,
    colorBuffer,
    pointBuffer,
    flagBuffer,
    camera,
    projectionTF
  ) {
    const { universe } = this.props;
    if (!this.reglCanvas || !universe) return;
    const cameraTF = camera.view();
    const projView = mat3.multiply(mat3.create(), projectionTF, cameraTF);
    const { width, height } = this.reglCanvas;
    regl.poll();
    regl.clear({
      depth: 1,
      color: [1, 1, 1, 1]
    });
    drawPoints({
      distance: camera.distance(),
      color: colorBuffer,
      position: pointBuffer,
      flag: flagBuffer,
      count: this.count,
      projView,
      nPoints: universe.nObs,
      minViewportDimension: Math.min(width || 800, height || 600)
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
      projectionTF
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
    const { responsive } = this.props;

    return (
      <div id="graphWrapper">
        <div
          style={{
            zIndex: -9999,
            position: "fixed",
            top: this.graphPaddingTop,
            right: globals.leftSidebarWidth
          }}
        >
          <div id="graphAttachPoint" />
          <div style={{ padding: 0, margin: 0 }}>
            <canvas
              width={responsive.width - this.graphPaddingRightLeft}
              height={responsive.height - this.graphPaddingTop}
              data-testid="layout-graph"
              ref={canvas => {
                this.reglCanvas = canvas;
              }}
              onMouseDown={this.handleCanvasEvent}
              onMouseUp={this.handleCanvasEvent}
              onMouseMove={this.handleCanvasEvent}
              onDoubleClick={this.handleCanvasEvent}
              onWheel={this.handleCanvasEvent}
            />
          </div>
        </div>
      </div>
    );
  }
}

export default Graph;
