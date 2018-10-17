// jshint esversion: 6
import React from "react";
import _ from "lodash";
import { connect } from "react-redux";
import * as d3 from "d3";
import { interpolateGreys } from "d3-scale-chromatic";
import HeatmapRow from "./heatmapRow";

/**********************************
***********************************
***********************************
        DiffExp Heatmap
***********************************
***********************************
**********************************/

@connect(state => ({
  differential: state.differential,
  world: state.controls.world
}))
class Heatmap extends React.Component {
  render() {
    const { world, differential } = this.props;
    if (!differential.diffExp) {
      return <p>Select cells & compute differential to see heatmap</p>;
    }

    // summarize the information for display.
    const topGenes = _.map(differential.diffExp, val => ({
      varIndex: val[0],
      geneName: world.varAnnotations[val[0]].name,
      avgDiff: val[1],
      set1AvgExp: val[4],
      set2AvgExp: val[5]
    }));

    // average expression extent
    const extent = d3.extent(
      _.concat(
        _.map(differential.diffExp, val => val[4]),
        _.map(differential.diffExp, val => val[5])
      )
    );

    const greyColorScale = d3
      .scaleSequential()
      .domain(extent)
      .interpolator(interpolateGreys);

    return (
      <div>
        <div
          style={{
            display: "flex",
            justifyContent: "flex-start",
            width: 400,
            fontWeight: 700
          }}
        >
          <p style={{ marginRight: 110 }}>Gene</p>
          <p style={{ marginRight: 25 }}>1</p>
          <p style={{ marginRight: 20 }}>2</p>
          <p>ave diff</p>
        </div>
        {topGenes.map(g => {
          const { geneName, avgDiff, set1AvgExp, set2AvgExp } = g;
          return (
            <HeatmapRow
              key={geneName}
              gene={geneName}
              greyColorScale={greyColorScale}
              aveDiff={avgDiff}
              set1exp={set1AvgExp}
              set2exp={set2AvgExp}
            />
          );
        })}
      </div>
    );
  }
}

export default Heatmap;
