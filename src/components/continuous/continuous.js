import React from 'react';
import _ from "lodash";
import cells from "../../../data/GBM_metadata.js";
import drawParallelCoordinates from "./drawParallelCoordinates";
import createContinuousRanges from "./createContinuousRanges";

const Continuous = () => {

  drawParallelCoordinates(cells)

  return (
    <div>
      <h3> Continuous Metadata </h3>
    </div>
  )
};

export default Continuous;
