import React from "react";
import * as d3 from "d3";
import { connect, shallowEqual } from "react-redux";
import { mat3, vec2 } from "gl-matrix";
import _regl from "regl";
import memoize from "memoize-one";
import Async from "react-async";
import { Button } from "@blueprintjs/core";

import setupSVGandBrushElements from "./setupSVGandBrush";
import _camera from "../../util/camera";
import _drawPoints from "./drawPointsRegl";
import {
  createColorTable,
  createColorQuery,
} from "../../util/stateManager/colorHelpers";
import * as globals from "../../globals";

import GraphOverlayLayer from "./overlays/graphOverlayLayer";
import CentroidLabels from "./overlays/centroidLabels";
import actions from "../../actions";
import renderThrottle from "../../util/renderThrottle";

import {
  flagBackground,
  flagSelected,
  flagHighlight,
} from "../../util/glHelpers";

/*
Simple 2D transforms control all point painting.  There are three:
  * model - convert from underlying per-point coordinate to a layout.
    Currently used to move from data to webgl coordinate system.
  * camera - apply a 2D camera transformation (pan, zoom)
  * projection - apply any transformation required for screen size and layout
*/
// eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
function createProjectionTF(viewportWidth: any, viewportHeight: any) {
  /*
  the projection transform accounts for the screen size & other layout
  */
  const fractionToUse = 0.95; // fraction of min dimension to use
  const topGutterSizePx = 32; // top gutter for tools
  const bottomGutterSizePx = 32; // bottom gutter for tools
  const heightMinusGutter =
    viewportHeight - topGutterSizePx - bottomGutterSizePx;
  const minDim = Math.min(viewportWidth, heightMinusGutter);
  const aspectScale = [
    (fractionToUse * minDim) / viewportWidth,
    (fractionToUse * minDim) / viewportHeight,
  ];
  const m = mat3.create();
  mat3.fromTranslation(m, [
    0,
    (bottomGutterSizePx - topGutterSizePx) / viewportHeight / aspectScale[1],
  ]);
  // @ts-expect-error ts-migrate(2345) FIXME: Argument of type 'number[]' is not assignable to p... Remove this comment to see the full error message
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

// eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
type GraphState = any;

// @ts-expect-error ts-migrate(1238) FIXME: Unable to resolve signature of class decorator whe... Remove this comment to see the full error message
@connect((state) => ({
  // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
  annoMatrix: (state as any).annoMatrix,
  // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
  crossfilter: (state as any).obsCrossfilter,
  // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
  selectionTool: (state as any).graphSelection.tool,
  // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
  currentSelection: (state as any).graphSelection.selection,
  // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
  layoutChoice: (state as any).layoutChoice,
  // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
  graphInteractionMode: (state as any).controls.graphInteractionMode,
  // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
  colors: (state as any).colors,
  // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
  pointDilation: (state as any).pointDilation,
  // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
  genesets: (state as any).genesets.genesets,
}))
// eslint-disable-next-line @typescript-eslint/ban-types --- FIXME: disabled temporarily on migrate to TS.
class Graph extends React.Component<{}, GraphState> {
  // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any -- - FIXME: disabled temporarily on migrate to TS.
  static createReglState(canvas: any) {
    /*
        Must be created for each canvas
        */
    // setup canvas, webgl draw function and camera
    const camera = _camera(canvas);
    const regl = _regl(canvas);
    const drawPoints = _drawPoints(regl);
    // preallocate webgl buffers
    // @ts-expect-error ts-migrate(2554) FIXME: Expected 1 arguments, but got 0.
    const pointBuffer = regl.buffer();
    // @ts-expect-error ts-migrate(2554) FIXME: Expected 1 arguments, but got 0.
    const colorBuffer = regl.buffer();
    // @ts-expect-error ts-migrate(2554) FIXME: Expected 1 arguments, but got 0.
    const flagBuffer = regl.buffer();
    return {
      camera,
      regl,
      drawPoints,
      pointBuffer,
      colorBuffer,
      flagBuffer,
    };
  }

  // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any -- - FIXME: disabled temporarily on migrate to TS.
  static watchAsync(props: any, prevProps: any) {
    return !shallowEqual(props.watchProps, prevProps.watchProps);
  }

  // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
  cachedAsyncProps: any;

  // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
  reglCanvas: any;

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
    (crossfilter, _flagSelected, _flagUnselected) => {
      const x = crossfilter.fillByIsSelected(
        new Float32Array(crossfilter.size()),
        _flagSelected,
        _flagUnselected
      );
      return x;
    }
  );

  computeHighlightFlags = memoize(
    (nObs, pointDilationData, pointDilationLabel) => {
      const flags = new Float32Array(nObs);
      if (pointDilationData) {
        for (let i = 0, len = flags.length; i < len; i += 1) {
          if (pointDilationData[i] === pointDilationLabel) {
            flags[i] = flagHighlight;
          }
        }
      }
      return flags;
    }
  );

  computeColorByFlags = memoize((nObs, colorByData) => {
    const flags = new Float32Array(nObs);
    if (colorByData) {
      for (let i = 0, len = flags.length; i < len; i += 1) {
        const val = colorByData[i];
        if (typeof val === "number" && !Number.isFinite(val)) {
          flags[i] = flagBackground;
        }
      }
    }
    return flags;
  });

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
      const nObs = crossfilter.size();
      const flags = new Float32Array(nObs);
      const selectedFlags = this.computeSelectedFlags(
        crossfilter,
        flagSelected,
        0
      );
      const highlightFlags = this.computeHighlightFlags(
        nObs,
        pointDilationData,
        pointDilationLabel
      );
      const colorByFlags = this.computeColorByFlags(nObs, colorByData);
      for (let i = 0; i < nObs; i += 1) {
        flags[i] = selectedFlags[i] + highlightFlags[i] + colorByFlags[i];
      }
      return flags;
    }
  );

  // eslint-disable-next-line @typescript-eslint/ban-types --- FIXME: disabled temporarily on migrate to TS.
  constructor(props: {}) {
    super(props);
    const viewport = this.getViewportDimensions();
    this.reglCanvas = null;
    this.cachedAsyncProps = null;
    const modelTF = createModelTF();
    this.state = {
      toolSVG: null,
      tool: null,
      container: null,
      viewport,
      // projection
      camera: null,
      modelTF,
      // @ts-expect-error ts-migrate(2345) FIXME: Argument of type '[]' is not assignable to paramet... Remove this comment to see the full error message
      modelInvTF: mat3.invert([], modelTF),
      projectionTF: createProjectionTF(viewport.width, viewport.height),
      // regl state
      regl: null,
      drawPoints: null,
      pointBuffer: null,
      colorBuffer: null,
      flagBuffer: null,
      // component rendering derived state - these must stay synchronized
      // with the reducer state they were generated from.
      layoutState: {
        layoutDf: null,
        layoutChoice: null,
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

  // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types --- FIXME: disabled temporarily on migrate to TS.
  componentDidMount() {
    window.addEventListener("resize", this.handleResize);
  }

  // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/ban-types --- FIXME: disabled temporarily on migrate to TS.
  componentDidUpdate(prevProps: {}, prevState: GraphState) {
    const {
      // @ts-expect-error ts-migrate(2339) FIXME: Property 'selectionTool' does not exist on type 'R... Remove this comment to see the full error message
      selectionTool,
      // @ts-expect-error ts-migrate(2339) FIXME: Property 'currentSelection' does not exist on type 'R... Remove this comment to see the full error message
      currentSelection,
      // @ts-expect-error ts-migrate(2339) FIXME: Property 'graphInteractionMode' does not exist on type 'R... Remove this comment to see the full error message
      graphInteractionMode,
    } = this.props;
    const { toolSVG, viewport } = this.state;
    const hasResized =
      prevState.viewport.height !== viewport.height ||
      prevState.viewport.width !== viewport.width;
    let stateChanges = {};
    if (
      (viewport.height && viewport.width && !toolSVG) || // first time init
      hasResized || //  window size has changed we want to recreate all SVGs
      // @ts-expect-error ts-migrate(2339) FIXME: Property 'selectionTool' does not exist on type '{... Remove this comment to see the full error message
      selectionTool !== prevProps.selectionTool || // change of selection tool
      // @ts-expect-error ts-migrate(2339) FIXME: Property 'graphInteractionMode' does not exist on ... Remove this comment to see the full error message
      prevProps.graphInteractionMode !== graphInteractionMode // lasso/zoom mode is switched
    ) {
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
      // @ts-expect-error ts-migrate(2339) FIXME: Property 'currentSelection' does not exist on type... Remove this comment to see the full error message
      currentSelection !== prevProps.currentSelection ||
      // @ts-expect-error ts-migrate(2339) FIXME: Property 'graphInteractionMode' does not exist on ... Remove this comment to see the full error message
      graphInteractionMode !== prevProps.graphInteractionMode ||
      // @ts-expect-error ts-migrate(2339) FIXME: Property 'toolSVG' does not exist on type '{}'.
      stateChanges.toolSVG
    ) {
      const { tool, container } = this.state;
      this.selectionToolUpdate(
        // @ts-expect-error ts-migrate(2339) FIXME: Property 'tool' does not exist on type '{}'.
        stateChanges.tool ? stateChanges.tool : tool,
        // @ts-expect-error ts-migrate(2339) FIXME: Property 'container' does not exist on type '{}'.
        stateChanges.container ? stateChanges.container : container
      );
    }
    if (Object.keys(stateChanges).length > 0) {
      // eslint-disable-next-line react/no-did-update-set-state --- Preventing update loop via stateChanges and diff checks
      this.setState(stateChanges);
    }
  }

  // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types --- FIXME: disabled temporarily on migrate to TS.
  componentWillUnmount() {
    window.removeEventListener("resize", this.handleResize);
  }

  // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types --- FIXME: disabled temporarily on migrate to TS.
  handleResize = () => {
    const { state } = this.state;
    const viewport = this.getViewportDimensions();
    const projectionTF = createProjectionTF(viewport.width, viewport.height);
    this.setState({
      ...state,
      viewport,
      projectionTF,
    });
  };

  // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any -- - FIXME: disabled temporarily on migrate to TS.
  handleCanvasEvent = (e: any) => {
    const { camera, projectionTF } = this.state;
    if (e.type !== "wheel") e.preventDefault();
    if (camera.handleEvent(e, projectionTF)) {
      this.renderCanvas();
      // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
      this.setState((state: any) => {
        return { ...state, updateOverlay: !state.updateOverlay };
      });
    }
  };

  // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types --- FIXME: disabled temporarily on migrate to TS.
  handleBrushDragAction() {
    /*
          event describing brush position:
          @-------|
          |       |
          |       |
          |-------@
        */
    // ignore programatically generated events
    // @ts-expect-error ts-migrate(2339) FIXME: Property 'event' does not exist on type 'typeof im... Remove this comment to see the full error message
    if (d3.event.sourceEvent === null || !d3.event.selection) return;
    // @ts-expect-error ts-migrate(2339) FIXME: Property 'dispatch' does not exist on type 'Readon... Remove this comment to see the full error message
    const { dispatch, layoutChoice } = this.props;
    // @ts-expect-error ts-migrate(2339) FIXME: Property 'event' does not exist on type 'typeof im... Remove this comment to see the full error message
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

  // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types --- FIXME: disabled temporarily on migrate to TS.
  handleBrushStartAction() {
    // Ignore programatically generated events.
    // @ts-expect-error ts-migrate(2339) FIXME: Property 'event' does not exist on type 'typeof im... Remove this comment to see the full error message
    if (!d3.event.sourceEvent) return;
    // @ts-expect-error ts-migrate(2339) FIXME: Property 'dispatch' does not exist on type 'Readon... Remove this comment to see the full error message
    const { dispatch } = this.props;
    dispatch(actions.graphBrushStartAction());
  }

  // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types --- FIXME: disabled temporarily on migrate to TS.
  handleBrushEndAction() {
    // Ignore programatically generated events.
    // @ts-expect-error ts-migrate(2339) FIXME: Property 'event' does not exist on type 'typeof im... Remove this comment to see the full error message
    if (!d3.event.sourceEvent) return;
    /*
        coordinates will be included if selection made, null
        if selection cleared.
        */
    // @ts-expect-error ts-migrate(2339) FIXME: Property 'dispatch' does not exist on type 'Readon... Remove this comment to see the full error message
    const { dispatch, layoutChoice } = this.props;
    // @ts-expect-error ts-migrate(2339) FIXME: Property 'event' does not exist on type 'typeof im... Remove this comment to see the full error message
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

  // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types --- FIXME: disabled temporarily on migrate to TS.
  handleBrushDeselectAction() {
    // @ts-expect-error ts-migrate(2339) FIXME: Property 'dispatch' does not exist on type 'Readon... Remove this comment to see the full error message
    const { dispatch, layoutChoice } = this.props;
    dispatch(actions.graphBrushDeselectAction(layoutChoice.current));
  }

  // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types --- FIXME: disabled temporarily on migrate to TS.
  handleLassoStart() {
    // @ts-expect-error ts-migrate(2339) FIXME: Property 'dispatch' does not exist on type 'Readon... Remove this comment to see the full error message
    const { dispatch, layoutChoice } = this.props;
    // @ts-expect-error ts-migrate(2554) FIXME: Expected 0 arguments, but got 1.
    dispatch(actions.graphLassoStartAction(layoutChoice.current));
  }

  // when a lasso is completed, filter to the points within the lasso polygon
  // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any -- - FIXME: disabled temporarily on migrate to TS.
  handleLassoEnd(polygon: any) {
    const minimumPolygonArea = 10;
    // @ts-expect-error ts-migrate(2339) FIXME: Property 'dispatch' does not exist on type 'Readon... Remove this comment to see the full error message
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
          // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
          polygon.map((xy: any) => this.mapScreenToPoint(xy))
        )
      );
    }
  }

  // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types --- FIXME: disabled temporarily on migrate to TS.
  handleLassoCancel() {
    // @ts-expect-error ts-migrate(2339) FIXME: Property 'dispatch' does not exist on type 'Readon... Remove this comment to see the full error message
    const { dispatch, layoutChoice } = this.props;
    dispatch(actions.graphLassoCancelAction(layoutChoice.current));
  }

  // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types --- FIXME: disabled temporarily on migrate to TS.
  handleLassoDeselectAction() {
    // @ts-expect-error ts-migrate(2339) FIXME: Property 'dispatch' does not exist on type 'Readon... Remove this comment to see the full error message
    const { dispatch, layoutChoice } = this.props;
    dispatch(actions.graphLassoDeselectAction(layoutChoice.current));
  }

  // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types --- FIXME: disabled temporarily on migrate to TS.
  handleDeselectAction() {
    // @ts-expect-error ts-migrate(2339) FIXME: Property 'selectionTool' does not exist on type 'R... Remove this comment to see the full error message
    const { selectionTool } = this.props;
    if (selectionTool === "brush") this.handleBrushDeselectAction();
    if (selectionTool === "lasso") this.handleLassoDeselectAction();
  }

  // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any -- - FIXME: disabled temporarily on migrate to TS.
  handleOpacityRangeChange(e: any) {
    // @ts-expect-error ts-migrate(2339) FIXME: Property 'dispatch' does not exist on type 'Readon... Remove this comment to see the full error message
    const { dispatch } = this.props;
    dispatch({
      type: "change opacity deselected cells in 2d graph background",
      data: e.target.value,
    });
  }

  // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any -- - FIXME: disabled temporarily on migrate to TS.
  setReglCanvas = (canvas: any) => {
    this.reglCanvas = canvas;
    this.setState({
      ...Graph.createReglState(canvas),
    });
  };

  // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types --- FIXME: disabled temporarily on migrate to TS.
  getViewportDimensions = () => {
    // @ts-expect-error ts-migrate(2339) FIXME: Property 'viewportRef' does not exist on type 'Rea... Remove this comment to see the full error message
    const { viewportRef } = this.props;
    return {
      height: viewportRef.clientHeight,
      width: viewportRef.clientWidth,
    };
  };

  // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types --- FIXME: disabled temporarily on migrate to TS.
  createToolSVG = () => {
    /*
        Called from componentDidUpdate. Create the tool SVG, and return any
        state changes that should be passed to setState().
        */
    // @ts-expect-error ts-migrate(2339) FIXME: Property 'selectionTool' does not exist on type 'R... Remove this comment to see the full error message
    const { selectionTool, graphInteractionMode } = this.props;
    const { viewport } = this.state;
    /* clear out whatever was on the div, even if nothing, but usually the brushes etc */
    const lasso = d3.select("#lasso-layer");
    if (lasso.empty()) return {}; // still initializing
    lasso.selectAll(".lasso-group").remove();
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

  // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any -- - FIXME: disabled temporarily on migrate to TS.
  fetchAsyncProps = async (props: any) => {
    const {
      annoMatrix,
      colors: colorsProp,
      layoutChoice,
      crossfilter,
      pointDilation,
      viewport,
    } = props.watchProps;
    const { modelTF } = this.state;
    const [layoutDf, colorDf, pointDilationDf] = await this.fetchData(
      annoMatrix,
      layoutChoice,
      colorsProp,
      pointDilation
    );
    const { currentDimNames } = layoutChoice;
    const X = layoutDf.col(currentDimNames[0]).asArray();
    const Y = layoutDf.col(currentDimNames[1]).asArray();
    const positions = this.computePointPositions(X, Y, modelTF);
    const colorTable = this.updateColorTable(colorsProp, colorDf);
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
    const { width, height } = viewport;
    return {
      positions,
      colors,
      flags,
      width,
      height,
    };
  };

  // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types --- FIXME: disabled temporarily on migrate to TS.
  async fetchData(
    // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any -- - FIXME: disabled temporarily on migrate to TS.
    annoMatrix: any,
    // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any -- - FIXME: disabled temporarily on migrate to TS.
    layoutChoice: any,
    // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any -- - FIXME: disabled temporarily on migrate to TS.
    colors: any,
    // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any -- - FIXME: disabled temporarily on migrate to TS.
    pointDilation: any
  ) {
    /*
        fetch all data needed.  Includes:
          - the color by dataframe
          - the layout dataframe
          - the point dilation dataframe
        */
    const { metadataField: pointDilationAccessor } = pointDilation;
    const promises = [];
    // layout
    promises.push(annoMatrix.fetch("emb", layoutChoice.current));
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

  // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any -- - FIXME: disabled temporarily on migrate to TS.
  brushToolUpdate(tool: any, container: any) {
    /*
        this is called from componentDidUpdate(), so be very careful using
        anything from this.state, which may be updated asynchronously.
        */
    // @ts-expect-error ts-migrate(2339) FIXME: Property 'currentSelection' does not exist on type... Remove this comment to see the full error message
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
                // @ts-expect-error ts-migrate(7053) FIXME: Element implicitly has an 'any' type because expre... Remove this comment to see the full error message
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

  // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any -- - FIXME: disabled temporarily on migrate to TS.
  lassoToolUpdate(tool: any) {
    /*
        this is called from componentDidUpdate(), so be very careful using
        anything from this.state, which may be updated asynchronously.
        */
    // @ts-expect-error ts-migrate(2339) FIXME: Property 'currentSelection' does not exist on type... Remove this comment to see the full error message
    const { currentSelection } = this.props;
    if (currentSelection.mode === "within-polygon") {
      /*
            if there is a current selection, make sure the lasso tool matches
            */
      // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
      const polygon = currentSelection.polygon.map((p: any) =>
        this.mapPointToScreen(p)
      );
      tool.move(polygon);
    } else {
      tool.reset();
    }
  }

  // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any -- - FIXME: disabled temporarily on migrate to TS.
  selectionToolUpdate(tool: any, container: any) {
    /*
        this is called from componentDidUpdate(), so be very careful using
        anything from this.state, which may be updated asynchronously.
        */
    // @ts-expect-error ts-migrate(2339) FIXME: Property 'selectionTool' does not exist on type 'R... Remove this comment to see the full error message
    const { selectionTool } = this.props;
    switch (selectionTool) {
      case "brush":
        this.brushToolUpdate(tool, container);
        break;
      case "lasso":
        // @ts-expect-error ts-migrate(2554) FIXME: Expected 1 arguments, but got 2.
        this.lassoToolUpdate(tool, container);
        break;
      default:
        /* punt? */
        break;
    }
  }

  // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any -- - FIXME: disabled temporarily on migrate to TS.
  mapScreenToPoint(pin: any) {
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

  // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any -- - FIXME: disabled temporarily on migrate to TS.
  mapPointToScreen(xyCell: any) {
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

  // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any -- - FIXME: disabled temporarily on migrate to TS.
  updateReglAndRender(asyncProps: any, prevAsyncProps: any) {
    const { positions, colors, flags, height, width } = asyncProps;
    this.cachedAsyncProps = asyncProps;
    const { pointBuffer, colorBuffer, flagBuffer } = this.state;
    let needToRenderCanvas = false;
    if (height !== prevAsyncProps?.height || width !== prevAsyncProps?.width) {
      needToRenderCanvas = true;
    }
    if (positions !== prevAsyncProps?.positions) {
      pointBuffer({ data: positions, dimension: 2 });
      needToRenderCanvas = true;
    }
    if (colors !== prevAsyncProps?.colors) {
      colorBuffer({ data: colors, dimension: 3 });
      needToRenderCanvas = true;
    }
    if (flags !== prevAsyncProps?.flags) {
      flagBuffer({ data: flags, dimension: 1 });
      needToRenderCanvas = true;
    }
    if (needToRenderCanvas) this.renderCanvas();
  }

  // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any -- - FIXME: disabled temporarily on migrate to TS.
  updateColorTable(colors: any, colorDf: any) {
    // @ts-expect-error ts-migrate(2339) FIXME: Property 'annoMatrix' does not exist on type 'Read... Remove this comment to see the full error message
    const { annoMatrix } = this.props;
    const { schema } = annoMatrix;
    /* update color table state */
    if (!colors || !colorDf) {
      return createColorTable(
        null, // default mode
        null,
        null,
        schema,
        null
      );
    }
    const { colorAccessor, userColors, colorMode } = colors;
    return createColorTable(
      colorMode,
      colorAccessor,
      colorDf,
      schema,
      userColors
    );
  }

  // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any -- - FIXME: disabled temporarily on migrate to TS.
  createColorByQuery(colors: any) {
    // @ts-expect-error ts-migrate(2339) FIXME: Property 'annoMatrix' does not exist on type 'Read... Remove this comment to see the full error message
    const { annoMatrix, genesets } = this.props;
    const { schema } = annoMatrix;
    const { colorMode, colorAccessor } = colors;
    return createColorQuery(colorMode, colorAccessor, schema, genesets);
  }

  // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types --- FIXME: disabled temporarily on migrate to TS.
  renderPoints(
    // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any -- - FIXME: disabled temporarily on migrate to TS.
    regl: any,
    // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any -- - FIXME: disabled temporarily on migrate to TS.
    drawPoints: any,
    // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any -- - FIXME: disabled temporarily on migrate to TS.
    colorBuffer: any,
    // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any -- - FIXME: disabled temporarily on migrate to TS.
    pointBuffer: any,
    // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any -- - FIXME: disabled temporarily on migrate to TS.
    flagBuffer: any,
    // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any -- - FIXME: disabled temporarily on migrate to TS.
    camera: any,
    // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any -- - FIXME: disabled temporarily on migrate to TS.
    projectionTF: any
  ) {
    // @ts-expect-error ts-migrate(2339) FIXME: Property 'annoMatrix' does not exist on type 'Read... Remove this comment to see the full error message
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

  // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types --- FIXME: disabled temporarily on migrate to TS.
  render() {
    const {
      // @ts-expect-error ts-migrate(2339) FIXME: Property 'graphInteractionMode' does not exist on ... Remove this comment to see the full error message
      graphInteractionMode,
      // @ts-expect-error ts-migrate(2339) FIXME: Property 'annoMatrix' does not exist on ... Remove this comment to see the full error message
      annoMatrix,
      // @ts-expect-error ts-migrate(2339) FIXME: Property 'colors' does not exist on ... Remove this comment to see the full error message
      colors,
      // @ts-expect-error ts-migrate(2339) FIXME: Property 'layoutChoice' does not exist on ... Remove this comment to see the full error message
      layoutChoice,
      // @ts-expect-error ts-migrate(2339) FIXME: Property 'pointDilation' does not exist on ... Remove this comment to see the full error message
      pointDilation,
      // @ts-expect-error ts-migrate(2339) FIXME: Property 'crossfilter' does not exist on ... Remove this comment to see the full error message
      crossfilter,
    } = this.props;
    const { modelTF, projectionTF, camera, viewport, regl } = this.state;
    const cameraTF = camera?.view()?.slice();
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
          // @ts-expect-error ts-migrate(2322) FIXME: Type '{ children: Element; width: any; height: any... Remove this comment to see the full error message
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
          ref={this.setReglCanvas}
          onMouseDown={this.handleCanvasEvent}
          onMouseUp={this.handleCanvasEvent}
          onMouseMove={this.handleCanvasEvent}
          onDoubleClick={this.handleCanvasEvent}
          onWheel={this.handleCanvasEvent}
        />

        <Async
          watchFn={Graph.watchAsync}
          promiseFn={this.fetchAsyncProps}
          watchProps={{
            annoMatrix,
            colors,
            layoutChoice,
            pointDilation,
            crossfilter,
            viewport,
          }}
        >
          <Async.Pending initial>
            {/* eslint-disable-next-line @typescript-eslint/no-use-before-define --- StillLoading used before defined */}
            <StillLoading
              displayName={layoutChoice.current}
              width={viewport.width}
              height={viewport.height}
            />
          </Async.Pending>
          <Async.Rejected>
            {(error) => (
              // eslint-disable-next-line @typescript-eslint/no-use-before-define --- ErrorLoading used before defined
              <ErrorLoading
                displayName={layoutChoice.current}
                error={error}
                width={viewport.width}
                height={viewport.height}
              />
            )}
          </Async.Rejected>
          <Async.Fulfilled>
            {(asyncProps) => {
              if (regl && !shallowEqual(asyncProps, this.cachedAsyncProps)) {
                this.updateReglAndRender(asyncProps, this.cachedAsyncProps);
              }
              return null;
            }}
          </Async.Fulfilled>
        </Async>
      </div>
    );
  }
}

// eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
const ErrorLoading = ({ displayName, error, width, height }: any) => {
  console.log(error); // log to console as this is an unepected error
  return (
    <div
      style={{
        position: "fixed",
        fontWeight: 500,
        top: height / 2,
        left: globals.leftSidebarWidth + width / 2 - 50,
      }}
    >
      <span>{`Failure loading ${displayName}`}</span>
    </div>
  );
};

// eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
const StillLoading = ({ displayName, width, height }: any) => {
  /*
  Render a busy/loading indicator
  */
  return (
    <div
      style={{
        position: "fixed",
        fontWeight: 500,
        top: height / 2,
        width,
      }}
    >
      <div
        style={{
          display: "flex",
          justifyContent: "center",
          justifyItems: "center",
          alignItems: "center",
        }}
      >
        <Button minimal loading intent="primary" />
        <span style={{ fontStyle: "italic" }}>Loading {displayName}</span>
      </div>
    </div>
  );
};

export default Graph;
