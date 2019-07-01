// jshint esversion: 6
import React from "react";
import { connect } from "react-redux";
import * as d3 from "d3";
import {
  Popover,
  PopoverInteractionKind,
  Position,
  Button,
  Classes
} from "@blueprintjs/core";

@connect()
class Occupancy extends React.Component {
  createHistogram = bins => {
    if (!this.canvas) return;
    const width = 100;
    const height = 11;

    const xScale = d3
      .scaleLinear()
      .domain([0, bins.length])
      .range([0, width]);

    const largestBin = Math.max(...bins);

    const yScale = d3
      .scaleLinear()
      .domain([0, largestBin])
      .range([0, height]);

    const ctx = this.canvas.getContext("2d");

    ctx.fillStyle = "#bbb";

    let x;
    let y;

    const rectWidth = width / bins.length;

    for (let i = 0, { length } = bins; i < length; i += 1) {
      x = xScale(i);
      y = yScale(bins[i]);
      ctx.fillRect(x, height - y, rectWidth, y);
    }
  };

  createOccupancyStack = stacks => {
    const height = 11;
    const ctx = this.canvas?.getContext("2d");

    if (!ctx) return;

    for (let i = 0, { length } = stacks; i < length; i += 1) {
      const { rectWidth, offset, fill } = stacks[i];
      ctx.fillStyle = fill;
      ctx.fillRect(offset, 0, rectWidth, height);
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

    this.canvas?.getContext("2d").clearRect(0, 0, width, height);

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
      <Popover
        interactionKind={PopoverInteractionKind.HOVER}
        hoverOpenDelay={1000}
        position={Position.LEFT}
        modifiers={{
          preventOverflow: { enabled: false },
          hide: { enabled: false }
        }}
        lazy
        usePortal
        popoverClassName={Classes.POPOVER_CONTENT_SIZING}
      >
        <canvas
          className="bp3-popover-targer"
          style={{
            marginRight: 5,
            width,
            height
          }}
          width={width}
          height={height}
          ref={ref => {
            this.canvas = ref;
            if (colorByIsCatagoricalData) this.createOccupancyStack(stacks);
            else this.createHistogram(occupancy);
          }}
        />
        <div key="text" style={{ fontFamily: "Roboto", fontSize: "14px" }}>
          <p>
            The x axis scale for all mini histograms are the same as the active
            color scale: <strong>{colorAccessor}</strong>
          </p>
          <p>
            The y axis scale for each mini histogram is normalized to its
            largest value
          </p>
        </div>
      </Popover>
    );
  }
}

export default Occupancy;
