import React from 'react';
import _ from "lodash";
import styles from "./graph.css";
import {setupGraphElements, drawGraph} from "./drawGraph";
import SectionHeader from "../framework/sectionHeader";
import { connect } from "react-redux";
import ColorControl from "../controls/color";

@connect((state) => {

  const vertices = state.cells.cells && state.cells.cells.data.graph ? state.cells.cells.data.graph : null;
  const ranges = state.cells.cells && state.cells.cells.data.ranges ? state.cells.cells.data.ranges : null;
  const metadata = state.cells.cells && state.cells.cells.data.metadata ? state.cells.cells.data.metadata : null;

  return {
    ranges,
    vertices,
    metadata,
    color: state.controls.color,
    continuousSelection: state.controls.continuousSelection,
  }
})
class Graph extends React.Component {

  constructor(props) {
    super(props);
    this.state = {
      drawn: false,
      svg: null,
      context: null,
    };
  }

  componentWillReceiveProps(nextProps) {
    /* maybe should do a check here to confirm ref exists and pass it? */
    if (
      this.state.context &&
      nextProps.vertices
      // nextProps.expressions &&
      // nextProps.expressionsCountsMap &&
    ) {
      drawGraph(
        nextProps.vertices,
        this.state.context,
        nextProps.expressionsCountsMap,
        nextProps.color,
        nextProps.ranges, /* assumption that this exists if vertices does both are on cells */
        nextProps.metadata,
        nextProps.continuousSelection, /* continuousSelected should probably inform a global 'is active' array rather than be consumed so specifically here */
      )
    }
  }

  componentDidMount() {
    const {svg, context} = setupGraphElements();
    this.setState({svg, context});
  }

  render() {
    return (
      <div id="graphWrapper">
        <SectionHeader text="Graph"/>
        <ColorControl/>
        <div id="graphAttachPoint"> </div>
      </div>
    )
  }
};

export default Graph;
