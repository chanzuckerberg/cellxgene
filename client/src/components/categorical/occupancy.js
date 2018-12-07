// jshint esversion: 6
import React from "react";
import _ from "lodash";
import { connect } from "react-redux";
import * as d3 from "d3";

@connect()
class Occupancy extends React.Component {
  render() {
    const {
      occupancy,
      colorScale,
      metadataField,
      colorAccessor,
      schema
    } = this.props;
    const width = 100;
    const height = 11;

    const categories = _.filter(schema.annotations.obs, {
      name: colorAccessor
    })[0].categories;

    const x = d3
      .scaleLinear()
      /* get all the keys d[1] as an array, then find the sum */
      .domain([0, d3.sum(Array.from(occupancy, d => d[1]))])
      .range([0, width]);

    const stack = [];
    let currentOffset = 0;
    occupancy.forEach((value, key) => {
      const scaledValue = x(value);
      stack.push({
        key,
        value,
        rectWidth: scaledValue,
        offset: currentOffset,
        fill: colorScale(categories.indexOf(key))
      });
      currentOffset += scaledValue;
    });

    return (
      <svg
        style={{
          marginRight: 5,
          width,
          height
        }}
      >
        {_.map(stack, d => (
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
