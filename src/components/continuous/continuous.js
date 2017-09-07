import React from 'react';
import _ from "lodash";
import cells from "../../../data/GBM_metadata.js";
import drawParallelCoordinates from "./drawParallelCoordinates";
import createContinuousRanges from "./createContinuousRanges";
import styles from './parallelCoordinates.css';

import {
  margin,
  width,
  height
} from "./util";

class Continuous extends React.Component {

  constructor(props) {
    super(props);
    this.state = {
      ctx: null,

    };
  }

  componentDidMount() {
    drawParallelCoordinates(cells)
  }

  render() {
    return (
      <div id="parcoords_wrapper" style={{marginTop: 50}}>
        <h3> Continuous Metadata </h3>
        <div style={{marginBottom: 30}}>
          Color By:
          <button style={{marginLeft: 10}}>Cluster_2d_color</button>
          <button style={{marginLeft: 10}}>Cluster_CNV_color</button>
          <button style={{marginLeft: 10}}>Location.color</button>
          <button style={{marginLeft: 10}}>Sample.name.color</button>
          <button style={{marginLeft: 10}}>Sample.type.color</button>
          <button style={{marginLeft: 10}}>Selection.color</button>
          <button style={{marginLeft: 10}}>housekeeping_cluster_color</button>
          <button style={{marginLeft: 10}}>recluster_myeloid</button>
          <button style={{marginLeft: 10}}>recluster_myeloid_color</button>
        </div>
        <div
          className={styles.parcoords}
          id="parcoords"
          style={{
            width:  width + margin.left + margin.right + "px",
            height: height + margin.top + margin.bottom + "px"
          }}></div>
        <pre id="parcoords_output" className={styles.pre}></pre>
      </div>
    )
  }
};

export default Continuous;
