import React from 'react';
import _ from "lodash";
import * as globals from "../../globals";
import styles from "./graph.css";
import {setupGraphElements, drawGraphUsingRenderQueue} from "./drawGraph";
import SectionHeader from "../framework/sectionHeader";
import { connect } from "react-redux";
import actions from "../../actions";
import drawScatterplotRegl from "./drawGraphRegl";

import mat4 from 'gl-mat4';
import fit from 'canvas-fit';
import _camera from '../../util/camera.js'
import _regl from 'regl'
import _drawPoints from './drawPointsRegl'

// import data generators
const generatePoints = function (count) {
  return Array(count).fill().map(function () {
    return [
      Math.random() - 0.5,
      Math.random() - 0.5
    ]
  })
}

const generateColors = function (count) {
  return Array(count).fill().map(function () {
    return [
      Math.random(),
      Math.random(),
      Math.random()
    ]
  })
}
// set constants
const count = 1993


@connect((state) => {

  const vertices = state.cells.cells && state.cells.cells.data.graph ? state.cells.cells.data.graph : null;
  const ranges = state.cells.cells && state.cells.cells.data.ranges ? state.cells.cells.data.ranges : null;
  const metadata = state.cells.cells && state.cells.cells.data.metadata ? state.cells.cells.data.metadata : null;

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
    opacityForDeselectedCells: state.controls.opacityForDeselectedCells,
  }
})
class Graph extends React.Component {

  constructor(props) {
    super(props);
    this.state = {
      drawn: false,
      svg: null,
      ctx: null,
      brush: null,
    };
  }
  componentDidMount() {
    // const {svg, ctx} = setupGraphElements(
    //   this.handleBrushSelectAction.bind(this),
    //   this.handleBrushDeselectAction.bind(this)
    // );
    // this.setState({svg, ctx});

    // setup canvas and camera
    const camera = _camera(this.reglCanvas, {scale: true, rotate: false});
    const regl = _regl(this.reglCanvas)

    const drawPoints = _drawPoints(regl)

    // preallocate buffers
    const pointBuffer = regl.buffer(generatePoints(count))
    const colorBuffer = regl.buffer(generateColors(count))


    regl.frame(({viewportWidth, viewportHeight}) => {

      regl.clear({
        depth: 1,
        color: [1, 1, 1, 1]
      })

      drawPoints({
        distance: camera.distance,
        color: colorBuffer,
        position: pointBuffer,
        count: count,
        view: camera.view(),
        scale: viewportHeight / viewportWidth
      })

      camera.tick()
    })

    this.setState({
      regl,
      pointBuffer,
      colorBuffer,
    })

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
    if (
      this.state.regl &&
      nextProps.vertices
    ) {

      const _currentCellSelectionMap = _.keyBy(nextProps.currentCellSelection, "CellName"); /* move me to the reducer */

      const positions = [];
      const colors = [];

      _.each(nextProps.currentCellSelection, (cell, i) => {
        if (nextProps.graphMap[cell["CellName"]]) { /* fails silently, sometimes this is undefined, in which case the graph array should be shorter than the cell array, check in reducer */
          positions.push([
            nextProps.graphMap[cell["CellName"]][0],
            nextProps.graphMap[cell["CellName"]][1]
          ])

          const scaleRGB = (input) => {
            const outputMax = 1;
            const outputMin = 0;

            const inputMax = 255;
            const inputMin = 0;

            const percent = (input - inputMin) / (inputMax - inputMin);
            return percent * (outputMax - outputMin) + outputMin;
          }

          let c = cell["__color__"];

          if (c[0] !== "#") {
            const _c = c.replace(/[^\d,.]/g, '').split(',');
            colors.push([
              scaleRGB(+_c[0]),
              scaleRGB(+_c[1]),
              scaleRGB(+_c[2])
            ])
          } else {
            var parsedHex = /^#?([a-f\d]{2})([a-f\d]{2})([a-f\d]{2})$/i.exec(c);
            colors.push([
              scaleRGB(parseInt(parsedHex[1], 16)),
              scaleRGB(parseInt(parsedHex[2], 16)),
              scaleRGB(parseInt(parsedHex[3], 16))
            ]);
          }

        }
      })

      this.state.pointBuffer(positions)
      this.state.colorBuffer(colors)
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
    }

    brushCoords.dx = brushCoords.southeastX - brushCoords.northwestX;
    brushCoords.dy = brushCoords.southeastY - brushCoords.northwestY;

    this.props.dispatch({
      type: "graph brush selection change",
      brushCoords
    })
  }
  handleBrushDeselectAction() {
    if (!d3.event.selection) {
      this.props.dispatch({
        type: "graph brush deselect"
      })
    }
  }
  handleOpacityRangeChange(e) {
    this.props.dispatch({
      type: "change opacity deselected cells in 2d graph background",
      data: e.target.value
    })
  }

  render() {
    return (
      <div
        id="graphWrapper"
        style={{
          height: 540, /* move this to globals */
          backgroundColor: "white",
          borderRadius: 3,
          boxShadow: "3px 4px 13px 0px rgba(201,201,201,1)",
        }}>
        <div style={{padding: 10}}>
          <button
            onClick={() => { this.props.dispatch(actions.regraph()) }}
            style={{
              fontSize: 12,
              fontWeight: 400,
              color: "white",
              padding: "10px 20px",
              backgroundColor: globals.darkGreen,
              border: "none",
              cursor: "pointer",
            }}
          >
          Regraph present selection
          </button>
          <span style={{ marginLeft: 20}}>
            deselected opacity
          </span>
          <input
            style={{position: "relative", top: 3, left: 10}}
            type="range"
            onChange={this.handleOpacityRangeChange.bind(this)}
            min={0}
            max={1}
            step="0.01"
            />
        </div>
        <div
          id="graphAttachPoint"
          >
        </div>
        <div style={{padding: 0, margin: 0}}>
          <canvas width={1500} height={800} style={{border: "1px solid black"}}  ref={(canvas) => { this.reglCanvas = canvas}}/>
        </div>
      </div>
    )
  }
};

export default Graph;
