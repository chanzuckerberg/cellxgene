import React from 'react';
import _ from "lodash";
// import styles from "./joy.css";
// import drawGraph from "./drawGraph";
// import joyParser from "./joyParser";

class Graph extends React.Component {

  constructor(props) {
    super(props);
    this.state = {

    };
  }

  componentWillReceiveProps(nextProps) {
    if (nextProps.data) {
      console.log('graph data 55', nextProps.data)
      // drawGraph(graphParser(nextProps.data));
    }
  }

  componentDidMount() {

  }

  render() {
    return (
      <div id="joyplot_wrapper" style={{marginTop: 50}}>
        <h3> Graph </h3>
        <div id="joyplot"> </div>
      </div>
    )
  }
};

export default Graph;
