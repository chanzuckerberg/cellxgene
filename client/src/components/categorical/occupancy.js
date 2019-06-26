// jshint esversion: 6
import React from "react";
import { connect } from "react-redux";
import * as d3 from "d3";

@connect()
class Occupancy extends React.Component {
  createKDE = occupancy => {
    if (!this.svg) return;
    const width = 100;
    const height = 11;

    const xScale = d3
      .scaleLinear()
      .domain([0, occupancy.length])
      .range([0, width]);

    const largestBin = occupancy.reduce((max, length) => Math.max(max, length));

    const yScale = d3
      .scaleLinear()
      .domain([0, largestBin])
      .range([0, height]);

    const ctx = this.svg.getContext("2d");

    ctx.fillStyle = "#bbb";

    let x;
    let y;

    const rectWidth = width / occupancy.length;

    for (let i = 0, { length } = occupancy; i < length; i += 1) {
      x = xScale(i);
      y = yScale(occupancy[i]);
      ctx.fillRect(x, height - y, rectWidth, y);
    }
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
    } else {
      const groupBy = world.obsAnnotations.col(metadataField);

      const col =
        world.obsAnnotations.col(colorAccessor) ||
        world.varData.col(colorAccessor);

      const range = col.summarize();

      occupancyMap = col.histogram(
        50,
        [range.min, range.max],
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
    }

    return (
      <canvas
        style={{
          marginRight: 5,
          width,
          height
        }}
        width={width}
        height={height}
        ref={ref => {
          this.svg = ref;
          if (!colorByIsCatagoricalData) this.createKDE(occupancy);
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
      </canvas>
    );
  }
}

export default Occupancy;
