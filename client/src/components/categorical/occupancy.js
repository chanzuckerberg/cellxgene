// jshint esversion: 6
import React from "react";
import { connect } from "react-redux";
import * as d3 from "d3";

@connect()
class Occupancy extends React.Component {
  render() {
    const {
      occupancy,
      colorScale,
      categoricalSelection,
      colorAccessor,
      schema
    } = this.props;
    const width = 100;
    const height = 11;

    const categories = schema.annotations.obsByName[colorAccessor]?.categories;

    const x = d3
      .scaleLinear()
      /* get all the keys d[1] as an array, then find the sum */
      .domain([0, d3.sum(Array.from(occupancy, d => d[1]))])
      .range([0, width]);

    let currentOffset = 0;

    const stacks = categoricalSelection[colorAccessor].categoryValues.map(d => {
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

    return (
      <svg
        style={{
          marginRight: 5,
          width,
          height
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
