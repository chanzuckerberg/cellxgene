import React from 'react';
import _ from "lodash";
import * as globals from "../../globals";
import styles from "./graph.css";
import {setupGraphElements, drawGraphUsingRenderQueue} from "./drawGraph";
import SectionHeader from "../framework/sectionHeader";
import { connect } from "react-redux";
import actions from "../../actions";

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
  componentWillReceiveProps(nextProps) {
    /* maybe should do a check here to confirm ref exists and pass it? */
    if (
      this.state.ctx &&
      nextProps.vertices
      // nextProps.expressions &&
      // nextProps.expressionsCountsMap &&
    ) {
      /* clear canvas */
      this.state.ctx.clearRect(0, 0, globals.graphWidth, globals.graphHeight);

      drawGraphUsingRenderQueue(
        this.state.ctx,
        nextProps.expressionsCountsMap,
        nextProps.colorAccessor,
        nextProps.ranges, /* assumption that this exists if vertices does both are on cells */
        nextProps.metadata,
        nextProps.currentCellSelection,
        nextProps.graphBrushSelection,
        nextProps.colorScale,
        nextProps.graphMap,
        nextProps.opacityForDeselectedCells,
      )
    }
  }

  componentDidMount() {
    const {svg, ctx} = setupGraphElements(
      this.handleBrushSelectAction.bind(this),
      this.handleBrushDeselectAction.bind(this)
    );
    this.setState({svg, ctx});
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
        <div style={{position: "relative", left: 20, top: 20}}>
          <button
            onClick={() => { this.props.dispatch(actions.regraph()) }}
            style={{
              fontSize: 14,
              fontWeight: 400,
              color: "#41633C",
              padding: "10px 20px",
              backgroundColor: "#B9D1B5",
              border: "none",
              borderRadius: 3,
              cursor: "pointer",
            }}
          >
          Regraph
          </button>
          <span style={{ marginLeft: 20}}>
            background opacity
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
        <div id="graphAttachPoint"> </div>
      </div>
    )
  }
};

export default Graph;
