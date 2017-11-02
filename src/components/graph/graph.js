import React from 'react';
import _ from "lodash";
import styles from "./graph.css";
import {setupGraphElements, drawGraph} from "./drawGraph";
import SectionHeader from "../framework/sectionHeader";
import { connect } from "react-redux";

@connect((state) => {

  const vertices = state.cells.cells && state.cells.cells.data.graph ? state.cells.cells.data.graph : null;

  return {
    vertices
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
    console.log('le grpah', nextProps)
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
        nextProps.expressionsCountsMap
      );
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
        <div id="graphAttachPoint"> </div>
      </div>
    )
  }
};

export default Graph;
