// jshint esversion: 6
import React from "react";
import { connect } from "react-redux";
import * as d3 from "d3";
import {
  Popover,
  PopoverInteractionKind,
  Position,
  Classes
} from "@blueprintjs/core";

@connect()
class Occupancy extends React.Component {
  _WIDTH = 100;

  _HEIGHT = 11;

  createHistogram = () => {
    /* 
      Knowing that colorScale is based off continous data, 
      createHistogram fetches the continous data in relation to the cells releveant to the catagory value.
      It then seperates that data into 50 bins for drawing the mini-histogram
    */
    const {
      world,
      metadataField,
      colorAccessor,
      category,
      categoryIndex
    } = this.props;

    if (!this.canvas) return;

    const groupBy = world.obsAnnotations.col(metadataField);

    const col =
      world.obsAnnotations.col(colorAccessor) ||
      world.varData.col(colorAccessor);

    const range = col.summarize();

    const histogramMap = col.histogram(
      50,
      [range.min, range.max],
      groupBy
    ); /* Because the signature changes we really need different names for histogram to differentiate signatures  */

    const bins = histogramMap.get(category.categoryValues[categoryIndex]);

    const xScale = d3
      .scaleLinear()
      .domain([0, bins.length])
      .range([0, this._WIDTH]);

    const largestBin = Math.max(...bins);

    const yScale = d3
      .scaleLinear()
      .domain([0, largestBin])
      .range([0, this._HEIGHT]);

    const ctx = this.canvas.getContext("2d");

    ctx.fillStyle = "#000";

    let x;
    let y;

    const rectWidth = this._WIDTH / bins.length;

    for (let i = 0, { length } = bins; i < length; i += 1) {
      x = xScale(i);
      y = yScale(bins[i]);
      ctx.fillRect(x, this._HEIGHT - y, rectWidth, y);
    }
  };

  createOccupancyStack = () => {
    /* 
      Knowing that the color scale is based off of catagorical data, 
      createOccupancyStack obtains a map showing the number if cells per colored value
      Using the colorScale a stack of colored bars is drawn representing the map
     */
    const {
      world,
      metadataField,
      colorAccessor,
      category,
      categoryIndex,
      schema,
      colorScale
    } = this.props;

    const ctx = this.canvas?.getContext("2d");

    if (!ctx) return;

    const groupBy = world.obsAnnotations.col(metadataField);
    const occupancyMap = world.obsAnnotations
      .col(colorAccessor)
      .histogram(groupBy);

    const occupancy = occupancyMap.get(category.categoryValues[categoryIndex]);

    const x = d3
      .scaleLinear()
      /* get all the keys d[1] as an array, then find the sum */
      .domain([0, d3.sum(Array.from(occupancy.values()))])
      .range([0, this._WIDTH]);
    const categories = schema.annotations.obsByName[colorAccessor]?.categories;

    let currentOffset = 0;
    const dfColumn = world.obsAnnotations.col(colorAccessor);
    const categoryValues = dfColumn.summarize().categories;

    let o;
    let scaledValue;
    let value;

    for (let i = 0, { length } = categoryValues; i < length; i += 1) {
      value = categoryValues[i];
      o = occupancy.get(value);
      scaledValue = x(o);
      ctx.fillStyle = o
        ? colorScale(categories.indexOf(value))
        : "rgb(255,255,255)";
      ctx.fillRect(currentOffset, 0, o ? scaledValue : 0, this._HEIGHT);
      currentOffset += o ? scaledValue : 0;
    }
  };

  render() {
    const { colorAccessor, categoricalSelection } = this.props;

    this.canvas?.getContext("2d").clearRect(0, 0, this._WIDTH, this._HEIGHT);

    const colorByIsCatagoricalData = !!categoricalSelection[colorAccessor];

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
        disabled={colorByIsCatagoricalData}
        popoverClassName={Classes.POPOVER_CONTENT_SIZING}
      >
        <canvas
          className="bp3-popover-targer"
          style={{
            marginRight: 5,
            width: this._WIDTH,
            height: this._HEIGHT,
            borderBottom: colorByIsCatagoricalData
              ? ""
              : "solid rgb(230, 230, 230) 0.25px"
          }}
          width={this._WIDTH}
          height={this._HEIGHT}
          ref={ref => {
            this.canvas = ref;
            if (colorByIsCatagoricalData) this.createOccupancyStack();
            else this.createHistogram();
          }}
        />
        <div key="text" style={{ fontFamily: "Roboto", fontSize: "14px" }}>
          <p style={{ margin: "0" }}>
            These histograms show the distribution of{" "}
            <strong>{colorAccessor}</strong> within each category. The x axis is
            the same for each histogram, while the y axis is scaled to the
            highest bin within each histogram.
          </p>
        </div>
      </Popover>
    );
  }
}

export default Occupancy;
