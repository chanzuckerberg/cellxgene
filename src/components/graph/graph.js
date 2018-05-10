// jshint esversion: 6
import React from "react";
import _ from "lodash";
import * as globals from "../../globals";
import styles from "./graph.css";
import { setupGraphElements } from "./drawGraph";
import SectionHeader from "../framework/sectionHeader";
import { connect } from "react-redux";
import actions from "../../actions";

import mat4 from "gl-mat4";
import fit from "canvas-fit";
import _camera from "../../util/camera.js";
import _regl from "regl";
import _drawPoints from "./drawPointsRegl";

import FaCrosshair from "react-icons/lib/fa/crosshairs";
import FaZoom from "react-icons/lib/fa/search-plus";
import FaSave from "react-icons/lib/fa/download";

/* https://bl.ocks.org/mbostock/9078690 - quadtree for onClick / hover selections */

@connect(state => {
  const vertices =
    state.cells.cells && state.cells.cells.data.graph
      ? state.cells.cells.data.graph
      : null;
  const ranges =
    state.cells.cells && state.cells.cells.data.ranges
      ? state.cells.cells.data.ranges
      : null;
  const metadata =
    state.cells.cells && state.cells.cells.data.metadata
      ? state.cells.cells.data.metadata
      : null;

  return {
    ranges,
    vertices,
    metadata,
    colorAccessor: state.controls.colorAccessor,
    colorScale: state.controls.colorScale,
    continuousSelection: state.controls.continuousSelection,
    graphMap: state.controls.graphMap,
    currentCellSelection: state.controls.currentCellSelection,
    graphBrushSelection: state.controls.graphBrushSelection,
    opacityForDeselectedCells: state.controls.opacityForDeselectedCells
  };
})
class Graph extends React.Component {
  constructor(props) {
    super(props);
    this.count = 0;
    this.inverse = mat4.identity([]);
    this.state = {
      drawn: false,
      svg: null,
      ctx: null,
      brush: null,
      mode: "brush"
    };
  }
  componentDidMount() {
    const { svg } = setupGraphElements(
      this.handleBrushSelectAction.bind(this),
      this.handleBrushDeselectAction.bind(this)
    );
    this.setState({ svg });

    // setup canvas and camera
    const camera = _camera(this.reglCanvas, { scale: true, rotate: false });
    const regl = _regl(this.reglCanvas);

    const drawPoints = _drawPoints(regl);

    // preallocate buffers
    const pointBuffer = regl.buffer();
    const colorBuffer = regl.buffer();
    const sizeBuffer = regl.buffer();

    regl.frame(({ viewportWidth, viewportHeight }) => {
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
        view: camera.view(),
        scale: viewportHeight / viewportWidth
      });

      var view = camera.view(); // get the camera matrix
      var projection = mat4.perspective(
        [],
        Math.PI / 2,
        context.viewportWidth * props.scale / context.viewportHeight,
        0.01,
        1000
      ); // get the projection matrix
      var combined = mat.multiply([], projection, view); // this is the matrix applied to the transform
      this.inverse = mat.invert([], combined); // this is the inverse

      camera.tick();
    });

    this.setState({
      regl,
      pointBuffer,
      colorBuffer,
      sizeBuffer
    });
  }
  componentWillReceiveProps(nextProps) {
    /* maybe should do a check here to confirm ref exists and pass it? */
    // if (
    //   this.state.ctx &&
    //   nextProps.vertices
    //   // nextProps.expressions &&
    //   // nextProps.expressionsCountsMap &&
    // ) {
    //   drawGraphUsingRenderQueue(
    //     this.state.ctx,
    //     nextProps.expressionsCountsMap,
    //     nextProps.colorAccessor,
    //     nextProps.ranges, /* assumption that this exists if vertices does both are on cells */
    //     nextProps.metadata,
    //     nextProps.currentCellSelection,
    //     nextProps.graphBrushSelection,
    //     nextProps.colorScale,
    //     nextProps.graphMap,
    //     nextProps.opacityForDeselectedCells,
    //   )
    // }
    if (this.state.regl && nextProps.vertices) {
      const _currentCellSelectionMap = _.keyBy(
        nextProps.currentCellSelection,
        "CellName"
      ); /* move me to the reducer */

      const positions = [];
      positions.length = nextProps.currentCellSelection.length;
      const colors = [];
      colors.length = nextProps.currentCellSelection.length;
      const sizes = [];
      sizes.length = nextProps.currentCellSelection.length;

      const glScaleX = d3
        .scaleLinear()
        .domain([0, 1])
        .range([-1, 1]); /* padding */

      const glScaleY = d3
        .scaleLinear()
        .domain([0, 1])
        .range([1, -1]); /* padding */

      /*
        Construct Vectors
      */
      _.each(nextProps.currentCellSelection, (cell, i) => {
        if (nextProps.graphMap[cell["CellName"]]) {
          /* fails silently, sometimes this is undefined, in which case the graph array should be shorter than the cell array, check in reducer */
          positions[i] = [
            glScaleX(nextProps.graphMap[cell["CellName"]][0]),
            glScaleY(nextProps.graphMap[cell["CellName"]][1])
          ];

          colors[i] = cell.__colorRGB__;
          sizes[i] = cell["__selected__"]
            ? 4
            : 0.2; /* make this a function of the number of total cells, including regraph */
        }
      });

      this.state.pointBuffer(positions);
      this.state.colorBuffer(colors);
      this.state.sizeBuffer(sizes);
      this.count = positions.length;
    }
  }
  handleBrushSelectAction() {
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
    const brushCoords = {
      northwestX: s[0][0],
      northwestY: s[0][1],
      southeastX: s[1][0],
      southeastY: s[1][1]
    };

    brushCoords.dx = brushCoords.southeastX - brushCoords.northwestX;
    brushCoords.dy = brushCoords.southeastY - brushCoords.northwestY;

    this.props.dispatch({
      type: "graph brush selection change",
      brushCoords
    });
  }
  handleBrushDeselectAction() {
    if (!d3.event.selection) {
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
      <div
        id="graphWrapper"
        style={{
          height: 1050 /* move this to globals */,
          backgroundColor: "white",
          borderRadius: 3,
          boxShadow: "3px 4px 13px 0px rgba(201,201,201,1)"
        }}
      >
        <div
          style={{
            padding: 10,
            display: "flex",
            justifyContent: "space-between",
            alignItems: "baseline"
          }}
        >
          <button
            onClick={() => {
              this.props.dispatch(actions.regraph());
            }}
            style={{
              fontSize: 12,
              fontWeight: 400,
              color: "white",
              padding: "10px 20px",
              backgroundColor: globals.brightBlue,
              border: "none",
              cursor: "pointer"
            }}
          >
            Regraph present selection
          </button>
          <div>
            <span style={{ marginRight: 10, fontSize: 12 }}>
              deselected opacity
            </span>
            <input
              style={{ position: "relative", top: 6, marginRight: 20 }}
              type="range"
              onChange={this.handleOpacityRangeChange.bind(this)}
              min={0}
              max={1}
              step="0.01"
            />
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
                fontSize: 12,
                fontWeight: 400,
                color: "white",
                padding: "10px 20px",
                backgroundColor: globals.brightBlue,
                border: "none",
                cursor: "pointer"
              }}
            >
              {" "}
              <FaSave style={{ display: "inline-block" }} /> csv url for present
              selection{" "}
            </button>
          </div>
        </div>
        <div
          style={{ display: this.state.mode === "brush" ? "inherit" : "none" }}
          id="graphAttachPoint"
        />
        <div style={{ padding: 0, margin: 0 }}>
          <canvas
            width={globals.graphWidth}
            height={globals.graphHeight}
            ref={canvas => {
              this.reglCanvas = canvas;
            }}
          />
        </div>
      </div>
    );
  }
}

export default Graph;
