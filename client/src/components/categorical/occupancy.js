// jshint esversion: 6
import React from "react";
import { connect } from "react-redux";
import * as d3 from "d3";

@connect()
class Occupancy extends React.Component {
  createKDE = occupancy => {
    const width = 100;
    const height = 11;
    const kernelDensityEstimator = (kernel, X) => {
      return V => {
        return X.map(x => {
          return [
            x,
            d3.mean(V, v => {
              return kernel(x - v);
            })
          ];
        });
      };
    };

    const kernelEpanechnikov = k => {
      return v => {
        return Math.abs((v /= k)) <= 1 ? (0.75 * (1 - v * v)) / k : 0;
      };
    };

    const x = d3
      .scaleLinear()
      .domain([0, occupancy.length])
      .range([0, width]);

    const y = d3
      .scaleLinear()
      .domain([-0.3585, 5.648] /* TODO: hardcoded values for Apod */)
      .range([0, height]);

    const density = kernelDensityEstimator(kernelEpanechnikov(1), x.ticks(50))(
      occupancy
    );
    const svg = d3.select(this.svg);

    console.log(density);

    svg
      .append("path")
      .datum(density)
      .attr("fill", "none")
      .attr("stroke", "#000")
      .attr("stroke-width", 1.5)
      .attr("stroke-linejoin", "round")
      .attr(
        "d",
        d3
          .line()
          .curve(d3.curveBasis)
          .x(d => {
            return x(d[0]);
          })
          .y(d => {
            return y(d[1]);
          })
      );
  };

  render() {
    const {
      colorScale,
      colorAccessor,
      schema,
      world,
      categoricalSelection,
      metadataField,
      categoryIndex,
      category
    } = this.props;
    const width = 100;
    const height = 11;

    let occupancyMap = null;

    const colorByIsCatagoricalData = !!categoricalSelection[colorAccessor];

    if (colorByIsCatagoricalData) {
      const groupBy = world.obsAnnotations.col(metadataField);
      occupancyMap = world.obsAnnotations.col(colorAccessor).histogram(groupBy);
    } else if (false /* is continuous obsannotation (n_counts) */) {
      // return;
    } else {
      const groupBy = world.obsAnnotations.col(metadataField);

      occupancyMap = world.varData
        .col(
          colorAccessor
        ) /* this is magic and col knows what kind of data is being used for histo. Could consider separate function names to make the fork explicit*/
        .histogram(
          50,
          [-0.3585, 5.648] /* TODO: hardcoded APOD range */,
          groupBy
        ); /* Because the signature changes we really need different names for histogram to differentiate signatures  */
    }

    const occupancy = occupancyMap.get(category.categoryValues[categoryIndex]);
    let stacks = [];

    if (colorByIsCatagoricalData) {
      const x = d3
        .scaleLinear()
        /* get all the keys d[1] as an array, then find the sum */
        .domain([0, d3.sum(Array.from(occupancy, d => d[1]))])
        .range([0, width]);
      const categories =
        schema.annotations.obsByName[colorAccessor]?.categories;

      let currentOffset = 0;
      const dfColumn = world.obsAnnotations.col(colorAccessor);
      const categoryValues = dfColumn.summarize().categories;
      stacks = categoryValues.map(d => {
        const o = occupancy.get(d);

        const scaledValue = x(o);

        const stackItem = {
          key: d,
          value: o || 0,
          rectWidth: o ? scaledValue : 0,
          offset: currentOffset,
          fill: o ? colorScale(categories.indexOf(d)) : "rgb(255,255,255)"
        };
        currentOffset += o ? scaledValue : 0;
        return stackItem;
      });
    } else {
    }

    return (
      <svg
        style={{
          marginRight: 5,
          width,
          height
        }}
        ref={ref => {
          this.svg = ref;
          this.createKDE(occupancy);
        }}
      >
        {stacks.map(d => (
          <rect
            key={d.key}
            width={d.rectWidth}
            height={height}
            x={d.offset}
            title={d.metadataField}
            fill={d.fill}
          />
        ))}
      </svg>
    );
  }
}

export default Occupancy;
