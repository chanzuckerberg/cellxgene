import React from "react";
import * as d3 from "d3";
import {
  Popover,
  PopoverInteractionKind,
  Position,
  Classes,
} from "@blueprintjs/core";

export default class MiniHistogram extends React.PureComponent {
  constructor(props) {
    super(props);
    this.canvasRef = React.createRef();
    this.state = {};
  }

  drawHistogram = () => {
    /* 
      Knowing that colorScale is based off continuous data, 
      createHistogram fetches the continuous data in relation to the cells relevant to the category value.
      It then separates that data into 50 bins for drawing the mini-histogram
    */
    const {
      world,
      metadataField,
      colorAccessor,
      categoryValue,
      width,
      height,
    } = this.props;
    const ctx = this.canvasRef.current.getContext("2d");

    ctx.clearRect(0, 0, width, height);
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

    const bins = histogramMap.has(categoryValue)
      ? histogramMap.get(categoryValue)
      : new Array(50).fill(0);

    const xScale = d3.scaleLinear().domain([0, bins.length]).range([0, width]);

    const largestBin = Math.max(...bins);

    const yScale = d3.scaleLinear().domain([0, largestBin]).range([0, height]);

    ctx.fillStyle = "#000";

    let x;
    let y;

    const rectWidth = width / bins.length;

    for (let i = 0, { length } = bins; i < length; i += 1) {
      x = xScale(i);
      y = yScale(bins[i]);
      ctx.fillRect(x, height - y, rectWidth, y);
    }
  };

  componentDidMount = () => {
    this.drawHistogram();
  };

  componentDidUpdate = (prevProps) => {
    const { colorAccessor } = this.props;
    if (prevProps.colorAccessor !== colorAccessor) this.drawHistogram();
  };

  render() {
    const { colorAccessor, categoryValue, width, height } = this.props;

    return (
      <Popover
        interactionKind={PopoverInteractionKind.HOVER_TARGET_ONLY}
        hoverOpenDelay={1500}
        hoverCloseDelay={200}
        position={Position.LEFT}
        modifiers={{
          preventOverflow: { enabled: false },
          hide: { enabled: false },
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
            height,
            borderBottom: "solid rgb(230, 230, 230) 0.25px",
          }}
          width={width}
          height={height}
          ref={this.canvasRef}
        />
        <div key="text" style={{ fontSize: "14px" }}>
          <p style={{ margin: "0" }}>
            This histograms shows the distribution of{" "}
            <strong>{colorAccessor}</strong> within{" "}
            <strong>{categoryValue}</strong>.
            <br />
            <br />
            The x axis is the same for each histogram, while the y axis is
            scaled to the largest bin within this histogram instead of the
            largest bin within the whole category.
          </p>
        </div>
      </Popover>
    );
  }
}
