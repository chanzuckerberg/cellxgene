// jshint esversion: 6
import React from "react";
import _ from "lodash";
import * as d3 from "d3";
import { connect } from "react-redux";
import mat4 from "gl-mat4";
import _regl from "regl";

import FaCrosshair from "react-icons/lib/fa/crosshairs";
import FaZoom from "react-icons/lib/fa/search-plus";

import * as globals from "../../globals";
import setupSVGandBrushElements from "./setupSVGandBrush";
import actions from "../../actions";
import _camera from "../../util/camera";
import _drawPoints from "./drawPointsRegl";
import { scaleLinear } from "../../util/scaleLinear";

/* https://bl.ocks.org/mbostock/9078690 - quadtree for onClick / hover selections */

@connect(state => ({
  world: state.controls.world,
  crossfilter: state.controls.crossfilter,
  responsive: state.responsive,
  colorRGB: _.get(state.controls, "colorRGB", null),
  opacityForDeselectedCells: state.controls.opacityForDeselectedCells,
  selectionUpdate: _.get(state.controls, "crossfilter.updateTime", null)
}))
class Graph extends React.Component {
  constructor(props) {
    super(props);
    this.count = 0;
    this.inverse = mat4.identity([]);
    this.graphPaddingTop = 100;
    this.graphPaddingBottom = 45;
    this.graphPaddingRight = 10;
    this.renderCache = {
      positions: null,
      colors: null
    };
    this.state = {
      svg: null,
      brush: null,
      mode: "brush"
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
      sizeBuffer
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

        for (
          let i = 0, { positions } = this.renderCache;
          i < cellCount;
          i += 1
        ) {
          positions[2 * i] = glScaleX(obsLayout.X[i]);
          positions[2 * i + 1] = glScaleY(obsLayout.Y[i]);
        }
        pointBuffer({
          data: this.renderCache.positions,
          dimension: 2
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
      (responsive.height && responsive.width && !this.state.svg)
    ) {
      /* clear out whatever was on the div, even if nothing, but usually the brushes etc */
      d3.select("#graphAttachPoint")
        .selectAll("svg")
        .remove();
      const { svg, brush, brushContainer } = setupSVGandBrushElements(
        this.handleBrushSelectAction.bind(this),
        this.handleBrushDeselectAction.bind(this),
        responsive,
        this.graphPaddingTop
      );
      this.setState({ svg, brush, brushContainer });
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

  handleBrushSelectAction() {
    /*
    This conditional handles procedural brush deselect. Brush emits
    an event on procedural deselect because it is move: null
    */

    const { camera } = this.state;
    const { dispatch, responsive } = this.props;

    if (d3.event.sourceEvent !== null) {
      /*
      No idea why d3 event scope works like this
      but apparently
      it does
      https://bl.ocks.org/EfratVil/0e542f5fc426065dd1d4b6daaa345a9f
    */
      const s = d3.event.selection;
      /*
      event describing brush position:
      @-------|
      |       |
      |       |
      |-------@
    */

      // compute inverse view matrix
      const inverse = mat4.invert([], camera.view());

      // transform screen coordinates -> cell coordinates
      const invert = pin => {
        const x = (2 * pin[0]) / (responsive.height - this.graphPaddingTop) - 1;
        const y =
          2 * (1 - pin[1] / (responsive.height - this.graphPaddingTop)) - 1;
        const pout = [
          x * inverse[14] + inverse[12],
          y * inverse[14] + inverse[13]
        ];
        return [(pout[0] + 1) / 2, (pout[1] + 1) / 2];
      };

      const brushCoords = {
        northwest: invert([s[0][0], s[0][1]]),
        southeast: invert([s[1][0], s[1][1]])
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

  handleOpacityRangeChange(e) {
    const { dispatch } = this.props;
    dispatch({
      type: "change opacity deselected cells in 2d graph background",
      data: e.target.value
    });
  }

  render() {
    const { dispatch, responsive } = this.props;
    const { mode } = this.state;
    return (
      <div id="graphWrapper">
        <div style={{ position: "fixed", right: 0, top: 0 }}>
          <div
            style={{
              padding: 10,
              display: "flex",
              justifyContent: "flex-end",
              alignItems: "baseline"
            }}
          >
            <button
              type="button"
              onClick={() => {
                dispatch(actions.resetGraph());
              }}
              style={{
                fontSize: 14,
                fontWeight: 400,
                color: "white",
                padding: "0px 10px",
                height: 30,
                marginRight: 10,
                borderRadius: 2,
                backgroundColor: globals.brightBlue,
                border: "none",
                cursor: "pointer"
              }}
            >
              reset graph
            </button>
            <button
              type="button"
              onClick={() => {
                dispatch(actions.regraph());
              }}
              style={{
                fontSize: 14,
                fontWeight: 400,
                color: "white",
                padding: "0px 10px",
                height: 30,
                marginRight: 10,
                borderRadius: 2,
                backgroundColor: globals.brightBlue,
                border: "none",
                cursor: "pointer"
              }}
            >
              regraph selection
            </button>
            <div>
              <span>
                <button
                  type="button"
                  onClick={() => {
                    this.setState({ mode: "brush" });
                  }}
                  style={{
                    cursor: "pointer",
                    border:
                      mode === "brush" ? "1px solid black" : "1px solid white",
                    backgroundColor: "white",
                    padding: 5,
                    marginRight: 10,
                    borderRadius: 3
                  }}
                >
                  {" "}
                  <FaCrosshair />{" "}
                </button>
                <button
                  type="button"
                  onClick={() => {
                    this.handleBrushDeselectAction();
                    this.restartReglLoop();
                    this.setState({ mode: "zoom" });
                  }}
                  style={{
                    cursor: "pointer",
                    border:
                      mode === "zoom" ? "1px solid black" : "1px solid white",
                    backgroundColor: "white",
                    padding: 5,
                    marginRight: 10,
                    borderRadius: 3
                  }}
                >
                  {" "}
                  <FaZoom />{" "}
                </button>
              </span>
            </div>
          </div>
        </div>
        <div
          style={{
            marginRight: 50,
            marginTop: 50,
            zIndex: -9999,
            position: "fixed",
            right: this.graphPaddingRight,
            bottom: this.graphPaddingBottom
          }}
        >
          <div
            style={{
              display: mode === "brush" ? "inherit" : "none"
            }}
            id="graphAttachPoint"
          />
          <div style={{ padding: 0, margin: 0 }}>
            <canvas
              width={responsive.height - this.graphPaddingTop}
              height={responsive.height - this.graphPaddingTop}
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
