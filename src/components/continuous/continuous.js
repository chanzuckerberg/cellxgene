import React from 'react';
import _ from "lodash";
import cells from "../../../data/GBM_metadata.js";
import drawParallelCoordinates from "./drawParallelCoordinates";
import styles from './parallelCoordinates.css';
import SectionHeader from "../framework/sectionHeader";

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
        <SectionHeader text="Continuous Metadata"/>
        <div
          className={styles.parcoords}
          id="parcoords"
          style={{
            width:  width + margin.left + margin.right + "px",
            height: height + margin.top + margin.bottom + "px"
          }}></div>
      </div>
    )
  }
};

export default Continuous;

/* text appended to parellel coords giving N instances of selection */
// <pre id="parcoords_output" className={styles.pre}></pre>

/*

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


*/
