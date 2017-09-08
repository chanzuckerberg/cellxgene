import React from 'react';
import _ from "lodash";
import styles from "./graph.css";
import {setupGraphElements, drawGraph} from "./drawGraph";

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
      nextProps.data &&
      this.state.context
    ) {
      drawGraph(
        nextProps.data.data,
        this.state.context
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
        <h3> Graph </h3>
        <div id="graphAttachPoint"> </div>
      </div>
    )
  }
};

export default Graph;
