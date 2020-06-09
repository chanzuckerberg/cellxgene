// jshint esversion: 6
import React from "react";
import * as d3 from "d3";

export default class MiniStackedBar extends React.PureComponent {
  constructor(props) {
    super(props);
    this.canvasRef = React.createRef();
  }

  drawStacks = () => {
    /* 
      Knowing that the color scale is based off of catagorical data, 
      createOccupancyStack obtains a map showing the number if cells per colored value
      Using the colorScale a stack of colored bars is drawn representing the map
     */
    const {
      world,
      metadataField,
      colorAccessor,
      categoryValue,
      colorScale,
      width,
      height,
    } = this.props;
    const { schema } = world;

    const ctx = this.canvasRef?.current.getContext("2d");

    ctx.clearRect(0, 0, width, height);

    const groupBy = world.obsAnnotations.col(metadataField);
    const occupancyMap = world.obsAnnotations
      .col(colorAccessor)
      .histogramCategorical(groupBy);

    const occupancy = occupancyMap.get(categoryValue);

    if (occupancy && occupancy.size > 0) {
      // not all categories have occupancy, so occupancy may be undefined.
      const x = d3
        .scaleLinear()
        /* get all the keys d[1] as an array, then find the sum */
        .domain([0, d3.sum(Array.from(occupancy.values()))])
        .range([0, width]);
      const categories =
        schema.annotations.obsByName[colorAccessor]?.categories;

      let currentOffset = 0;
      const dfColumn = world.obsAnnotations.col(colorAccessor);
      const categoryValues = dfColumn.summarizeCategorical().categories;

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
        ctx.fillRect(currentOffset, 0, o ? scaledValue : 0, height);
        currentOffset += o ? scaledValue : 0;
      }
    }
  };

  componentDidUpdate = (prevProps) => {
    const { colorAccessor } = this.props;
    if (colorAccessor !== prevProps.colorAccessor) this.drawStacks();
  };

  componentDidMount = () => {
    this.drawStacks();
  };

  render() {
    const { width, height } = this.props;
    const { canvas } = this;
    if (canvas) canvas.getContext("2d").clearRect(0, 0, width, height);

    return (
      <canvas
        className="bp3-popover-targer"
        style={{
          marginRight: 5,
          width,
          height,
        }}
        width={width}
        height={height}
        ref={this.canvasRef}
      />
    );
  }
}
