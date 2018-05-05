// jshint esversion: 6
import React from "react";
import _ from "lodash";
import styles from "./joy.css";
import drawJoy from "./drawJoy";
import joyParser from "./joyParser";

class Joy extends React.Component {
  constructor(props) {
    super(props);
    this.state = {};
  }

  componentWillReceiveProps(nextProps) {
    if (nextProps.data) {
      console.log("joyplot data 44", nextProps.data);
      drawJoy(joyParser(nextProps.data));
    }
  }

  componentDidMount() {}

  render() {
    return (
      <div id="joyplot_wrapper" style={{ marginTop: 50 }}>
        <h3> Joy </h3>
        <p>
          {" "}
          Cell expression distribution per gene & if differential expression,
          Ie., cells for cluster 5, top genes expressed by cluster 8
        </p>
        <div id="joyplot"> </div>
      </div>
    );
  }
}

export default Joy;
