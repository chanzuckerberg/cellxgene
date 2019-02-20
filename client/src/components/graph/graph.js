// jshint esversion: 6
import React from "react";
import _ from "lodash";
import * as d3 from "d3";
import { connect } from "react-redux";
import mat4 from "gl-mat4";
import _regl from "regl";
import { Button, AnchorButton, Tooltip } from "@blueprintjs/core";
import * as globals from "../../globals";
import setupSVGandBrushElements from "./setupSVGandBrush";
import actions from "../../actions";
import _camera from "../../util/camera";
import _drawPoints from "./drawPointsRegl";
import scaleLinear from "../../util/scaleLinear";

/* https://bl.ocks.org/mbostock/9078690 - quadtree for onClick / hover selections */

@connect(state => ({
  world: state.controls.world,
  universe: state.controls.universe,
  crossfilter: state.controls.crossfilter,
  responsive: state.responsive,
  colorRGB: _.get(state.controls, "colorRGB", null),
  opacityForDeselectedCells: state.controls.opacityForDeselectedCells,
  selectionUpdate: _.get(state.controls, "crossfilter.updateTime", null),
  resettingInterface: state.controls.resettingInterface
}))
class Graph extends React.Component {
  constructor(props) {
    super(props);
    this.count = 0;
    this.inverse = mat4.identity([]);
    this.graphPaddingTop = 0;
    this.graphPaddingBottom = 45;
    this.graphPaddingRight = globals.leftSidebarWidth;
    this.renderCache = {
      positions: null,
      colors: null
    };
    this.state = {
      svg: null,
      brush: null,
      mode: "lasso"
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
      reglRender
    });
  }

  componentDidUpdate(prevProps) {
    const {
      world,
      crossfilter,
      selectionUpdate,
      colorRGB,
      responsive
    } = this.props;
    const {
      reglRender,
      mode,
      regl,
      drawPoints,
      camera,
      pointBuffer,
      colorBuffer,
      sizeBuffer,
      svg
    } = this.state;

    if (reglRender && this.reglRenderState === "rendering" && mode !== "zoom") {
      reglRender.cancel();
      this.reglRenderState = "paused";
    }

    if (regl && world) {
      /* update the regl state */
      const { obsLayout } = world;
      const cellCount = crossfilter.size();

      // X/Y positions for each point - a cached value that only
      // changes if we have loaded entirely new cell data
      //
      if (
        !this.renderCache.positions ||
        selectionUpdate !== prevProps.selectionUpdate
      ) {
        if (!this.renderCache.positions) {
          this.renderCache.positions = new Float32Array(2 * cellCount);
        }

        const glScaleX = scaleLinear([0, 1], [-1, 1]);
        const glScaleY = scaleLinear([0, 1], [1, -1]);

        const offset = [d3.mean(obsLayout.X) - 0.5, d3.mean(obsLayout.Y) - 0.5];

        for (
          let i = 0, { positions } = this.renderCache;
          i < cellCount;
          i += 1
        ) {
          positions[2 * i] = glScaleX(obsLayout.X[i] - offset[0]);
          positions[2 * i + 1] = glScaleY(obsLayout.Y[i] - offset[1]);
        }
        pointBuffer({
          data: this.renderCache.positions,
          dimension: 2
        });

        this.setState({
          offset
        });
      }

      // Colors for each point - a cached value that only changes when
      // the cell metadata changes (done by updateCellColors middleware).
      // NOTE: this is a slightly pessimistic assumption, as the metadata
      // could have changed for some other reason, but for now color is
      // the only metadata that changes client-side.  If this is problematic,
      // we could add some sort of color-specific indicator to the app state.
      if (!this.renderCache.colors || colorRGB !== prevProps.colorRGB) {
        const rgb = colorRGB;
        if (!this.renderCache.colors) {
          this.renderCache.colors = new Float32Array(3 * rgb.length);
        }
        for (let i = 0, { colors } = this.renderCache; i < rgb.length; i += 1) {
          colors.set(rgb[i], 3 * i);
        }
        colorBuffer({ data: this.renderCache.colors, dimension: 3 });
      }

      // Sizes for each point - this is presumed to change each time the
      // component receives new props.   Almost always a true assumption, as
      // most property upates are due to changes driving a crossfilter
      // selection set change.
      //
      if (!this.renderCache.sizes) {
        this.renderCache.sizes = new Float32Array(cellCount);
      }

      crossfilter.fillByIsFiltered(this.renderCache.sizes, 4, 0.2);
      sizeBuffer({ data: this.renderCache.sizes, dimension: 1 });

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

    if (
      prevProps.responsive.height !== responsive.height ||
      prevProps.responsive.width !== responsive.width ||
      /* first time */
      (responsive.height && responsive.width && !svg)
    ) {
      /* clear out whatever was on the div, even if nothing, but usually the brushes etc */
      d3.select("#graphAttachPoint")
        .selectAll("svg")
        .remove();
      const { svg: newSvg, brush } = setupSVGandBrushElements(
        this.handleBrushSelectAction.bind(this),
        this.handleBrushDeselectAction.bind(this),
        responsive,
        this.graphPaddingRight,
        this.handleLassoStart.bind(this),
        this.handleLassoEnd.bind(this)
      );
      this.setState({ svg: newSvg, brush });
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

  invertPoint(pin) {
    const { responsive } = this.props;
    const { regl, camera, offset } = this.state;

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
      y * inverse[14] + inverse[13]
    ];
    return [(pout[0] + 1) / 2 + offset[0], (pout[1] + 1) / 2 + offset[1]];
  }

  handleBrushSelectAction() {
    /*
      This conditional handles procedural brush deselect. Brush emits
      an event on procedural deselect because it is move: null
    */
    /*
      event describing brush position:
      @-------|
      |       |
      |       |
      |-------@
    */
    /*
      No idea why d3 event scope works like this
      but apparently
      it does
      https://bl.ocks.org/EfratVil/0e542f5fc426065dd1d4b6daaa345a9f
    */
    const { dispatch } = this.props;

    if (d3.event.sourceEvent !== null) {
      const s = d3.event.selection;

      const brushCoords = {
        northwest: this.invertPoint([s[0][0], s[0][1]]),
        southeast: this.invertPoint([s[1][0], s[1][1]])
      };

      dispatch({
        type: "graph brush selection change",
        brushCoords
      });
    }
  }

  handleBrushDeselectAction() {
    const { dispatch } = this.props;
    const { svg, brush } = this.state;

    if (d3.event && !d3.event.selection) {
      dispatch({
        type: "graph brush deselect"
      });
    }

    if (!d3.event) {
      /*
      this line clears the brush procedurally, ie., zoom button clicked,
      not a click away from brush on svg
      */
      svg.select(".graph_brush").call(brush.move, null);
      dispatch({
        type: "graph brush deselect"
      });
    }
  }

  handleLassoStart() {
    const { dispatch } = this.props;
    // reset selected points when starting a new polygon
    // making it easier for the user to make the next selection
    dispatch({
      type: "lasso started"
    });
  }

  // when a lasso is completed, filter to the points within the lasso polygon
  handleLassoEnd(polygon) {
    const { dispatch } = this.props;

    dispatch({
      type: "lasso selection",
      polygon: polygon.map(xy => this.invertPoint(xy)) // transform the polygon
    });
  }

  handleOpacityRangeChange(e) {
    const { dispatch } = this.props;
    dispatch({
      type: "change opacity deselected cells in 2d graph background",
      data: e.target.value
    });
  }

  resetInterface() {
    const { dispatch } = this.props;
    dispatch({
      type: "interface reset started"
    });
    dispatch(actions.resetInterface());
  }

  render() {
    const {
      dispatch,
      responsive,
      crossfilter,
      resettingInterface
    } = this.props;

    const { mode } = this.state;
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
                disabled={
                  crossfilter &&
                  (crossfilter.countFiltered() === 0 ||
                    crossfilter.countFiltered() === crossfilter.size())
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
                disabled={
                  false
                  /* world && universe ? worldEqUniverse(world, universe) : false */
                }
                type="button"
                loading={resettingInterface}
                intent="warning"
                style={{ marginRight: 10 }}
                onClick={this.resetInterface.bind(this)}
              >
                reset
              </AnchorButton>
            </Tooltip>
            <div>
              <div className="bp3-button-group">
                <Tooltip content="Lasso selection" position="left">
                  <Button
                    type="button"
                    className="bp3-button bp3-icon-polygon-filter"
                    active={mode === "lasso"}
                    onClick={() => {
                      this.handleBrushDeselectAction();
                      // this.restartReglLoop();
                      this.setState({ mode: "lasso" });
                    }}
                    style={{
                      cursor: "pointer"
                    }}
                  />
                </Tooltip>
                <Tooltip content="Pan and zoom" position="left">
                  <Button
                    type="button"
                    className="bp3-button bp3-icon-zoom-in"
                    active={mode === "zoom"}
                    onClick={() => {
                      this.handleBrushDeselectAction();
                      this.restartReglLoop();
                      this.setState({ mode: "zoom" });
                    }}
                    style={{
                      cursor: "pointer"
                    }}
                  />
                </Tooltip>
              </div>
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
              display: mode === "lasso" ? "inherit" : "none"
            }}
            id="graphAttachPoint"
          />
          <div style={{ padding: 0, margin: 0 }}>
            <canvas
              width={responsive.width - this.graphPaddingRight}
              height={responsive.height - this.graphPaddingTop}
              data-testid="layout"
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
