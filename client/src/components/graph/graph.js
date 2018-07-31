// jshint esversion: 6
import React from "react";
import _ from "lodash";
import * as globals from "../../globals";
import styles from "./graph.css";
import { setupSVGandBrushElements } from "./setupSVGandBrush";
import SectionHeader from "../framework/sectionHeader";
import { connect } from "react-redux";
import actions from "../../actions";

import mat4 from "gl-mat4";
import fit from "canvas-fit";
import _camera from "../../util/camera.js";
import _regl from "regl";
import _drawPoints from "./drawPointsRegl";
import { scaleLinear } from "../../util/scaleLinear";

import FaCrosshair from "react-icons/lib/fa/crosshairs";
import FaZoom from "react-icons/lib/fa/search-plus";
import FaSave from "react-icons/lib/fa/download";

/* https://bl.ocks.org/mbostock/9078690 - quadtree for onClick / hover selections */

@connect(state => {
  return {
    cellsMetadata: state.controls.cellsMetadata,
    opacityForDeselectedCells: state.controls.opacityForDeselectedCells,
    responsive: state.responsive,
    crossfilter: state.controls.crossfilter
  };
})
class Graph extends React.Component {
  constructor(props) {
    super(props);
    this.count = 0;
    this.inverse = mat4.identity([]);
    this.graphPaddingTop = 100;
    this.renderCache = {
      positions: null,
      colors: null
    };
    this.state = {
      drawn: false,
      svg: null,
      ctx: null,
      brush: null,
      mode: "brush"
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
  restartReglLoop() {
    const reglRender = this.state.regl.frame(() => {
      this.reglDraw(
        this.state.regl,
        this.state.drawPoints,
        this.state.sizeBuffer,
        this.state.colorBuffer,
        this.state.pointBuffer,
        this.state.camera
      );
      this.state.camera.tick();
    });

    this.reglRenderState = "rendering";

    this.setState({
      reglRender
    });
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
  // static getDerivedStateFromProps(props, state) {
  //   console.log("getDerivedStateFromProps in graph.js");
  //   console.log("props", props);
  //   console.log("state", state);
  //   // console.log("this.props", this.props);
  // }
  componentWillReceiveProps(nextProps) {
    if (this.state.regl && nextProps.crossfilter) {
      /* update the regl state */
      const crossfilter = nextProps.crossfilter.cells;
      const cells = crossfilter.all();
      const cellCount = cells.length;

      // X/Y positions for each point - a cached value that only
      // changes if we have loaded entirely new cell data
      //
      if (
        !this.renderCache.positions ||
        this.props.crossfilter.cells != nextProps.crossfilter.cells
      ) {
        if (!this.renderCache.positions)
          this.renderCache.positions = new Float32Array(2 * cellCount);

        // d3.scaleLinear().domain([0,1]).range([-1,1])
        const glScaleX = scaleLinear([0, 1], [-1, 1]);
        // d3.scaleLinear().domain([0,1]).range([1,-1])
        const glScaleY = scaleLinear([0, 1], [1, -1]);

        for (
          let i = 0, positions = this.renderCache.positions;
          i < cellCount;
          i++
        ) {
          positions[2 * i] = glScaleX(cells[i].__x__);
          positions[2 * i + 1] = glScaleY(cells[i].__y__);
        }
        this.state.pointBuffer({
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
      if (
        !this.renderCache.colors ||
        this.props.cellsMetadata != nextProps.cellsMetadata
      ) {
        if (!this.renderCache.colors)
          this.renderCache.colors = new Float32Array(3 * cellCount);
        for (let i = 0, colors = this.renderCache.colors; i < cellCount; i++) {
          colors.set(cells[i].__colorRGB__, 3 * i);
        }
        this.state.colorBuffer({ data: this.renderCache.colors, dimension: 3 });
      }

      // Sizes for each point - this is presumed to change each time the
      // component receives new props.   Almost always a true assumption, as
      // most property upates are due to changes driving a crossfilter
      // selection set change.
      //
      if (
        !this.renderCache.sizes ||
        this.props.crossfilter.cells != nextProps.crossfilter.cells
      ) {
        this.renderCache.sizes = new Float32Array(cellCount);
      }
      crossfilter.fillByIsFiltered(this.renderCache.sizes, 4, 0.2);
      this.state.sizeBuffer({ data: this.renderCache.sizes, dimension: 1 });

      this.count = cellCount;

      this.state.regl._refresh();
      this.reglDraw(
        this.state.regl,
        this.state.drawPoints,
        this.state.sizeBuffer,
        this.state.colorBuffer,
        this.state.pointBuffer,
        this.state.camera
      );
    }

    if (
      /* invisibly handles the initial null vs integer case as well as resize events */
      nextProps.responsive.height !== this.props.responsive.height ||
      nextProps.responsive.width !== this.props.responsive.width
    ) {
      /* clear out whatever was on the div, even if nothing, but usually the brushes etc */
      d3.select("#graphAttachPoint")
        .selectAll("svg")
        .remove();
      const { svg, brush, brushContainer } = setupSVGandBrushElements(
        this.handleBrushSelectAction.bind(this),
        this.handleBrushDeselectAction.bind(this),
        nextProps.responsive,
        this.graphPaddingTop
      );
      this.setState({ svg, brush, brushContainer });
    }
  }
  componentDidUpdate() {
    if (
      this.state.reglRender &&
      this.reglRenderState === "rendering" &&
      this.state.mode !== "zoom"
    ) {
      this.state.reglRender.cancel();
      this.reglRenderState = "paused";
    }
  }
  handleBrushSelectAction() {
    /* This conditional handles procedural brush deselect. Brush emits an event on procedural deselect because it is move: null */
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
      const inverse = mat4.invert([], this.state.camera.view());

      // transform screen coordinates -> cell coordinates
      const invert = pin => {
        const x =
          (2 * pin[0]) / (this.props.responsive.height - this.graphPaddingTop) -
          1;
        const y =
          2 *
            (1 -
              pin[1] / (this.props.responsive.height - this.graphPaddingTop)) -
          1;
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

      this.props.dispatch({
        type: "graph brush selection change",
        brushCoords
      });
    }
  }
  handleBrushDeselectAction() {
    if (d3.event && !d3.event.selection) {
      this.props.dispatch({
        type: "graph brush deselect"
      });
    }

    if (!d3.event) {
      /* this line clears the brush procedurally, ie., zoom button clicked, not a click away from brush on svg */
      this.state.svg.select(".graph_brush").call(this.state.brush.move, null);
      this.props.dispatch({
        type: "graph brush deselect"
      });
    }
  }
  handleOpacityRangeChange(e) {
    this.props.dispatch({
      type: "change opacity deselected cells in 2d graph background",
      data: e.target.value
    });
  }

  render() {
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
              onClick={() => {
                this.props.dispatch(actions.resetGraph());
              }}
              style={{
                fontSize: 14,
                fontWeight: 700,
                color: "white",
                padding: "10px 20px",
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
              onClick={() => {
                this.props.dispatch(actions.regraph());
              }}
              style={{
                fontSize: 14,
                fontWeight: 700,
                color: "white",
                padding: "10px 20px",
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
              <span style={{ position: "relative", top: 3 }}>
                <button
                  onClick={() => {
                    this.setState({ mode: "brush" });
                  }}
                  style={{
                    cursor: "pointer",
                    border:
                      this.state.mode === "brush"
                        ? "1px solid black"
                        : "1px solid white",
                    backgroundColor: "white",
                    padding: 5,
                    borderRadius: 3
                  }}
                >
                  {" "}
                  <FaCrosshair />{" "}
                </button>
                <button
                  onClick={() => {
                    this.handleBrushDeselectAction();
                    this.restartReglLoop();
                    this.setState({ mode: "zoom" });
                  }}
                  style={{
                    cursor: "pointer",
                    border:
                      this.state.mode === "zoom"
                        ? "1px solid black"
                        : "1px solid white",
                    backgroundColor: "white",
                    padding: 5,
                    borderRadius: 3
                  }}
                >
                  {" "}
                  <FaZoom />{" "}
                </button>
              </span>
            </div>
            <div>
              <button
                style={{
                  fontSize: 14,
                  fontWeight: 700,
                  color: "white",
                  padding: "10px 20px",
                  backgroundColor: globals.lightGrey,
                  border: "none",
                  cursor: "pointer"
                }}
              >
                export
              </button>
            </div>
          </div>
        </div>
        <div
          style={{
            marginRight: 50,
            marginTop: 50,
            zIndex: -9999,
            position: "fixed",
            right: 20,
            bottom: 20
          }}
        >
          <div
            style={{
              display: this.state.mode === "brush" ? "inherit" : "none"
            }}
            id="graphAttachPoint"
          />
          <div style={{ padding: 0, margin: 0 }}>
            <canvas
              width={this.props.responsive.height - this.graphPaddingTop}
              height={this.props.responsive.height - this.graphPaddingTop}
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

// <span style={{ marginRight: 10, fontSize: 12 }}>
//   deselected opacity
// </span>
// <input
//   style={{ position: "relative", top: 6, marginRight: 20 }}
//   type="range"
//   onChange={this.handleOpacityRangeChange.bind(this)}
//   min={0}
//   max={1}
//   step="0.01"
// />
