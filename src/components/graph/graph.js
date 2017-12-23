import React from 'react';
import _ from "lodash";
import styles from "./graph.css";
import {setupGraphElements, drawGraph} from "./drawGraph";
import SectionHeader from "../framework/sectionHeader";
import { connect } from "react-redux";

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
    graphBrushSelection: state.controls.graphBrushSelection
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
      drawGraph(
        this.state.ctx,
        nextProps.expressionsCountsMap,
        nextProps.colorAccessor,
        nextProps.ranges, /* assumption that this exists if vertices does both are on cells */
        nextProps.metadata,
        nextProps.currentCellSelection,
        nextProps.graphBrushSelection,
        nextProps.colorScale,
        nextProps.graphMap,
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

  render() {
    return (
      <div id="graphWrapper" style={{height: 600 /* move this to globals */}}>
        <div id="graphAttachPoint"> </div>
      </div>
    )
  }
};

export default Graph;

// <SectionHeader text="Graph"/>
