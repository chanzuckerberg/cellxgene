import React from 'react';
import _ from "lodash";
import styles from "./graph.css";
import drawGraph from "./drawGraph";

class Graph extends React.Component {

  constructor(props) {
    super(props);
    this.state = {
      drawn: false,
    };
  }

  componentWillReceiveProps(nextProps) {
    /* maybe should do a check here to confirm ref exists and pass it? */
    if (nextProps.data && !this.state.drawn) {
      console.log('About to draw graph. Data:', nextProps.data)
      drawGraph(nextProps.data.data)
      this.setState({drawn: true});
    }
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
