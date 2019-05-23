// jshint esversion: 6
import React from "react";
import * as d3 from "d3";
import { connect } from "react-redux";
import mat4 from "gl-mat4";
import _regl from "regl";
import memoize from "memoize-one";
import {
  Button,
  AnchorButton,
  Tooltip,
  Popover,
  Menu,
  MenuItem,
  Position,
  NumericInput,
  Icon,
  RadioGroup,
  Radio
} from "@blueprintjs/core";

import * as globals from "../../globals";
import setupSVGandBrushElements from "./setupSVGandBrush";
import actions from "../../actions";
import _camera from "../../util/camera";
import _drawPoints from "./drawPointsRegl";
import scaleLinear from "../../util/scaleLinear";
import { World } from "../../util/stateManager";

/* https://bl.ocks.org/mbostock/9078690 - quadtree for onClick / hover selections */

@connect(state => ({
  world: state.world,
  universe: state.universe,
  crossfilter: state.crossfilter,
  clipPercentileMin: Math.round(100 * (state.world?.clipQuantiles?.min ?? 0)),
  clipPercentileMax: Math.round(100 * (state.world?.clipQuantiles?.max ?? 1)),
  responsive: state.responsive,
  colorRGB: state.colors.rgb,
  opacityForDeselectedCells: state.controls.opacityForDeselectedCells,
  resettingInterface: state.controls.resettingInterface,
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
  layoutChoice: state.layoutChoice
}))
class Graph extends React.Component {
  static isValidDigitKeyEvent(e) {
    /*
    Return true if this event is necessary to enter a percent number input.
    Return false if not.

    Returns true for events with keys: backspace, control, alt, meta, [0-9],
    or events that don't have a key.
    */
    if (e.key === null) return true;
    if (e.ctrlKey || e.altKey || e.metaKey) return true;

    // concept borrowed from blueprint's numericInputUtils:
    // keys that print a single character when pressed have a `key` name of
    // length 1. every other key has a longer `key` name (e.g. "Backspace",
    // "ArrowUp", "Shift"). since none of those keys can print a character
    // to the field--and since they may have important native behaviors
    // beyond printing a character--we don't want to disable their effects.
    const isSingleCharKey = e.key.length === 1;
    if (!isSingleCharKey) return true;

    const key = e.key.charCodeAt(0) - 48; /* "0" */
    return key >= 0 && key <= 9;
  }

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
    this.graphPaddingTop = 0;
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
      container: null,
      mode: "select",
      pendingClipPercentiles: null
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
    const transform = {
      glScaleX: scaleLinear([0, 1], [-1, 1]),
      glScaleY: scaleLinear([0, 1], [1, -1])
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

  componentDidUpdate(prevProps, prevState) {
    const { renderCache } = this;
    const {
      world,
      crossfilter,
      colorRGB,
      responsive,
      selectionTool,
      currentSelection,
      layoutChoice
    } = this.props;
    const { reglRender, mode, regl, svg } = this.state;
    let stateChanges = {};

    if (reglRender && this.reglRenderState === "rendering" && mode !== "zoom") {
      reglRender.cancel();
      this.reglRenderState = "paused";
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
      mode !== prevState.mode ||
      stateChanges.svg
    ) {
      const { tool, container, transform } = this.state;
      this.selectionToolUpdate(
        stateChanges.tool ? stateChanges.tool : tool,
        stateChanges.container ? stateChanges.container : container,
        stateChanges.transform ? stateChanges.transform : transform
      );
    }

    if (Object.keys(stateChanges).length > 0) {
      this.setState(stateChanges);
    }
  }

  isResetDisabled = () => {
    /*
    Reset should be disabled when all of the following are true:
      * nothing is selected in the crossfilter
      * world EQ universe
      * nothing is colored by
      * there are no userDefinedGenes or diffexpGenes displayed
      * scatterplot is not displayed
      * nothing in cellset1 or cellset2
      * clip percentiles are [0,100]
    */
    const {
      crossfilter,
      world,
      universe,
      userDefinedGenes,
      diffexpGenes,
      colorAccessor,
      scatterplotXXaccessor,
      scatterplotYYaccessor,
      celllist1,
      celllist2,
      clipPercentileMin,
      clipPercentileMax
    } = this.props;

    if (!crossfilter || !world || !universe) {
      return false;
    }
    const nothingSelected = crossfilter.countSelected() === crossfilter.size();
    const nothingColoredBy = !colorAccessor;
    const noGenes = userDefinedGenes.length === 0 && diffexpGenes.length === 0;
    const scatterNotDpl = !scatterplotXXaccessor || !scatterplotYYaccessor;
    const nothingInCellsets = !celllist1 && !celllist2;

    return (
      nothingSelected &&
      World.worldEqUniverse(world, universe) &&
      nothingColoredBy &&
      noGenes &&
      scatterNotDpl &&
      nothingInCellsets &&
      clipPercentileMax === 100 &&
      clipPercentileMin === 0
    );
  };

  resetInterface = () => {
    const { dispatch } = this.props;
    dispatch({
      type: "interface reset started"
    });
    dispatch(actions.resetInterface());
  };

  isClipDisabled = () => {
    /*
    return true if clip button should be disabled.
    */
    const { pendingClipPercentiles } = this.state;
    const clipPercentileMin = pendingClipPercentiles?.clipPercentileMin;
    const clipPercentileMax = pendingClipPercentiles?.clipPercentileMax;

    const { world } = this.props;
    const currentClipMin = 100 * world?.clipQuantiles?.min;
    const currentClipMax = 100 * world?.clipQuantiles?.max;

    // if you change this test, be careful with logic around
    // comparisons between undefined / NaN handling.
    const isDisabled =
      !(clipPercentileMin < clipPercentileMax) ||
      (clipPercentileMin === currentClipMin &&
        clipPercentileMax === currentClipMax);

    return isDisabled;
  };

  handleClipOnKeyPress = e => {
    /*
    allow only numbers, plus other critical keys which
    may be required to make a number
    */
    if (!Graph.isValidDigitKeyEvent(e)) {
      e.preventDefault();
    }
  };

  handleClipPercentileMinValueChange = v => {
    /*
    Ignore anything that isn't a legit number
    */
    if (!Number.isFinite(v)) return;

    const { pendingClipPercentiles } = this.state;
    const clipPercentileMax = pendingClipPercentiles?.clipPercentileMax;

    /*
    clamp to [0, currentClipPercentileMax]
    */
    if (v <= 0) v = 0;
    if (v > 100) v = 100;
    const clipPercentileMin = Math.round(v); // paranoia
    this.setState({
      pendingClipPercentiles: { clipPercentileMin, clipPercentileMax }
    });
  };

  handleClipPercentileMaxValueChange = v => {
    /*
    Ignore anything that isn't a legit number
    */
    if (!Number.isFinite(v)) return;

    const { pendingClipPercentiles } = this.state;
    const clipPercentileMin = pendingClipPercentiles?.clipPercentileMin;

    /*
    clamp to [0, 100]
    */
    if (v < 0) v = 0;
    if (v > 100) v = 100;
    const clipPercentileMax = Math.round(v); // paranoia

    this.setState({
      pendingClipPercentiles: { clipPercentileMin, clipPercentileMax }
    });
  };

  handleClipCommit = () => {
    const { dispatch } = this.props;
    const { pendingClipPercentiles } = this.state;
    const { clipPercentileMin, clipPercentileMax } = pendingClipPercentiles;
    const min = clipPercentileMin / 100;
    const max = clipPercentileMax / 100;
    dispatch({
      type: "set clip quantiles",
      clipQuantiles: { min, max }
    });
  };

  handleClipOpening = () => {
    const { clipPercentileMin, clipPercentileMax } = this.props;
    this.setState({
      pendingClipPercentiles: { clipPercentileMin, clipPercentileMax }
    });
  };

  handleClipClosing = () => {
    this.setState({ pendingClipPercentiles: null });
  };

  handleLayoutChoiceChange = e => {
    const { dispatch } = this.props;
    dispatch({
      type: "set layout choice",
      layoutChoice: e.currentTarget.value
    });
  };

  brushToolUpdate(tool, container, transform) {
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
          this.mapPointToScreen(
            currentSelection.brushCoords.northwest,
            transform
          ),
          this.mapPointToScreen(
            currentSelection.brushCoords.southeast,
            transform
          )
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

  lassoToolUpdate(tool, container, transform) {
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
        this.mapPointToScreen(p, transform)
      );
      tool.move(polygon);
    } else {
      tool.reset();
    }
  }

  selectionToolUpdate(tool, container, transform) {
    /*
    this is called from componentDidUpdate(), so be very careful using
    anything from this.state, which may be updated asynchronously.
    */
    const { selectionTool } = this.props;
    switch (selectionTool) {
      case "brush":
        this.brushToolUpdate(tool, container, transform);
        break;
      case "lasso":
        this.lassoToolUpdate(tool, container, transform);
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

    // compute inverse view matrix
    const inverse = mat4.invert([], camera.view());

    // transform screen coordinates -> cell coordinates
    const x = (2 * pin[0]) / (responsive.width - this.graphPaddingRight) - 1;
    const y = 2 * (1 - pin[1] / (responsive.height - this.graphPaddingTop)) - 1;
    const pout = [
      x * inverse[14] * aspect + inverse[12],
      -(y * inverse[14] + inverse[13])
    ];

    return [glScaleX.invert(pout[0]), glScaleY.invert(pout[1])];
  }

  mapPointToScreen(xyCell, transform) {
    /*
    Map an XY coordinate from cell/point domain to screen range.  Inverse
    of mapScreenToPoint()
    */
    const { responsive } = this.props;
    const { regl, camera } = this.state;
    const { glScaleX, glScaleY } = transform;

    const gl = regl._gl;

    // get aspect ratio
    const aspect = gl.drawingBufferWidth / gl.drawingBufferHeight;

    // compute inverse view matrix
    const inverse = mat4.invert([], camera.view());

    // variable names are choosen to reflect inverse of those used
    // in mapScreenToPoint().
    const pout = [glScaleX(xyCell[0]), glScaleY(xyCell[1])];
    const x = (pout[0] - inverse[12]) / aspect / inverse[14];
    const y = (-pout[1] - inverse[13]) / inverse[14];

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
    const {
      dispatch,
      responsive,
      crossfilter,
      resettingInterface,
      libraryVersions,
      undoDisabled,
      redoDisabled,
      selectionTool,
      clipPercentileMin,
      clipPercentileMax,
      layoutChoice
    } = this.props;
    const { mode, pendingClipPercentiles } = this.state;

    const clipMin =
      pendingClipPercentiles?.clipPercentileMin ?? clipPercentileMin;
    const clipMax =
      pendingClipPercentiles?.clipPercentileMax ?? clipPercentileMax;
    const activeClipClass =
      clipPercentileMin > 0 || clipPercentileMax < 100
        ? " bp3-intent-warning"
        : "";

    // constants used to create selection tool button
    let selectionTooltip;
    let selectionButtonClass;
    if (selectionTool === "brush") {
      selectionTooltip = "Brush selection";
      selectionButtonClass = "bp3-icon-select";
    } else {
      selectionTooltip = "Lasso selection";
      selectionButtonClass = "bp3-icon-polygon-filter";
    }

    return (
      <div id="graphWrapper">
        <div
          style={{
            position: "fixed",
            right: 0,
            top: 0
          }}
        >
          <div
            style={{
              padding: 10,
              display: "flex",
              justifyContent: "flex-end",
              alignItems: "baseline"
            }}
          >
            <Tooltip
              content="Show only metadata and cells which are currently selected"
              position="left"
            >
              <AnchorButton
                type="button"
                data-testid="subset-button"
                disabled={
                  crossfilter &&
                  (crossfilter.countSelected() === 0 ||
                    crossfilter.countSelected() === crossfilter.size())
                }
                style={{ marginRight: 10 }}
                onClick={() => {
                  dispatch(actions.regraph());
                  dispatch({ type: "increment graph render counter" });
                }}
              >
                subset to current selection
              </AnchorButton>
            </Tooltip>
            <Tooltip
              content="Reset cellxgene, clearing all selections"
              position="left"
            >
              <AnchorButton
                disabled={this.isResetDisabled()}
                type="button"
                loading={resettingInterface}
                intent="warning"
                style={{ marginRight: 10 }}
                onClick={this.resetInterface}
                data-testid="reset"
                data-testclass={`resetting-${resettingInterface}`}
              >
                reset
              </AnchorButton>
            </Tooltip>
            <div className="bp3-button-group">
              <Tooltip content={selectionTooltip} position="left">
                <Button
                  type="button"
                  data-testid="mode-lasso"
                  className={`bp3-button ${selectionButtonClass}`}
                  active={mode === "select"}
                  onClick={() => {
                    this.setState({ mode: "select" });
                  }}
                  style={{
                    cursor: "pointer"
                  }}
                />
              </Tooltip>
              <Tooltip content="Pan and zoom" position="left">
                <Button
                  type="button"
                  data-testid="mode-pan-zoom"
                  className="bp3-button bp3-icon-zoom-in"
                  active={mode === "zoom"}
                  onClick={() => {
                    this.restartReglLoop();
                    this.setState({ mode: "zoom" });
                  }}
                  style={{
                    cursor: "pointer"
                  }}
                />
              </Tooltip>
            </div>
            <div
              className="bp3-button-group"
              style={{
                marginLeft: 10
              }}
            >
              <Tooltip content="Undo" position="left">
                <AnchorButton
                  type="button"
                  className="bp3-button bp3-icon-undo"
                  disabled={undoDisabled}
                  onClick={() => {
                    dispatch({ type: "@@undoable/undo" });
                  }}
                  style={{
                    cursor: "pointer"
                  }}
                />
              </Tooltip>
              <Tooltip content="Redo" position="left">
                <AnchorButton
                  type="button"
                  className="bp3-button bp3-icon-redo"
                  disabled={redoDisabled}
                  onClick={() => {
                    dispatch({ type: "@@undoable/redo" });
                  }}
                  style={{
                    cursor: "pointer"
                  }}
                />
              </Tooltip>
            </div>

            <div
              className="bp3-button-group"
              style={{
                marginLeft: 10
              }}
            >
              <Popover
                target={
                  <Button
                    type="button"
                    data-testid="layout-choice"
                    className="bp3-button bp3-icon-heatmap"
                    style={{
                      cursor: "pointer"
                    }}
                  />
                }
                position={Position.BOTTOM_RIGHT}
                content={
                  <div
                    style={{
                      display: "flex",
                      justifyContent: "flex-start",
                      alignItems: "flex-start",
                      flexDirection: "column",
                      padding: 10
                    }}
                  >
                    <RadioGroup
                      label="Layout Choice"
                      onChange={this.handleLayoutChoiceChange}
                      selectedValue={layoutChoice.current}
                    >
                      {layoutChoice.available.map(name => (
                        <Radio label={name} value={name} key={name} />
                      ))}
                    </RadioGroup>
                  </div>
                }
              />
            </div>

            <div
              className="bp3-button-group"
              style={{
                marginLeft: 10
              }}
            >
              <Popover
                target={
                  <Button
                    type="button"
                    data-testid="visualization-settings"
                    className={`bp3-button bp3-icon-timeline-bar-chart ${activeClipClass}`}
                    style={{
                      cursor: "pointer"
                    }}
                  />
                }
                position={Position.BOTTOM_RIGHT}
                onOpening={this.handleClipOpening}
                onClosing={this.handleClipClosing}
                content={
                  <div
                    style={{
                      display: "flex",
                      justifyContent: "flex-start",
                      alignItems: "flex-start",
                      flexDirection: "column",
                      padding: 10
                    }}
                  >
                    <div>Clip all continuous values to percentile range</div>
                    <div
                      style={{
                        display: "flex",
                        justifyContent: "space-between",
                        alignItems: "center",
                        paddingTop: 5,
                        paddingBottom: 5
                      }}
                    >
                      <NumericInput
                        style={{ width: 50 }}
                        data-testid={"clip-min-input"}
                        onValueChange={this.handleClipPercentileMinValueChange}
                        onKeyPress={this.handleClipOnKeyPress}
                        value={clipMin}
                        min={0}
                        max={100}
                        fill={false}
                        minorStepSize={null}
                        rightElement={
                          <div style={{ padding: "4px 2px" }}>
                            <Icon
                              icon="percentage"
                              intent="primary"
                              iconSize={14}
                            />
                          </div>
                        }
                      />
                      <span style={{ marginRight: 5, marginLeft: 5 }}> - </span>
                      <NumericInput
                        style={{ width: 50 }}
                        data-testid={"clip-max-input"}
                        onValueChange={this.handleClipPercentileMaxValueChange}
                        onKeyPress={this.handleClipOnKeyPress}
                        value={clipMax}
                        min={0}
                        max={100}
                        fill={false}
                        minorStepSize={null}
                        rightElement={
                          <div style={{ padding: "4px 2px" }}>
                            <Icon
                              icon="percentage"
                              intent="primary"
                              iconSize={14}
                            />
                          </div>
                        }
                      />
                      <Button
                        type="button"
                        data-testid="clip-commit"
                        className="bp3-button"
                        disabled={this.isClipDisabled()}
                        style={{
                          cursor: "pointer",
                          marginRight: 5,
                          marginLeft: 5
                        }}
                        onClick={this.handleClipCommit}
                      >
                        Clip
                      </Button>
                    </div>
                  </div>
                }
              />
            </div>

            <div style={{ marginLeft: 10 }} className="bp3-button-group">
              <Popover
                content={
                  <Menu>
                    <MenuItem
                      href="https://chanzuckerberg.github.io/cellxgene/faq.html"
                      target="_blank"
                      icon="help"
                      text="FAQ"
                    />
                    <MenuItem
                      href="https://join-cellxgene-users.herokuapp.com/"
                      target="_blank"
                      icon="chat"
                      text="Chat"
                    />
                    <MenuItem
                      href="https://chanzuckerberg.github.io/cellxgene/"
                      target="_blank"
                      icon="book"
                      text="Docs"
                    />
                    <MenuItem
                      href="https://github.com/chanzuckerberg/cellxgene"
                      target="_blank"
                      icon="git-branch"
                      text="Github"
                    />
                    <MenuItem
                      target="_blank"
                      text={`cellxgene v${
                        libraryVersions && libraryVersions.cellxgene
                          ? libraryVersions.cellxgene
                          : null
                      }`}
                    />
                    <MenuItem text="MIT License" />
                  </Menu>
                }
                position={Position.BOTTOM_RIGHT}
              >
                <Button
                  type="button"
                  className="bp3-button bp3-icon-info-sign"
                  style={{
                    cursor: "pointer"
                  }}
                />
              </Popover>
            </div>
          </div>
        </div>
        <div
          style={{
            zIndex: -9999,
            position: "fixed",
            top: 0,
            right: 0
          }}
        >
          <div
            style={{
              display: mode === "select" ? "inherit" : "none"
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
