import React from "react";
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
  }

  drawHistogram = () => {
    const { xScale, yScale, bins, width, height } = this.props;

    if (!bins) return;

    const ctx = this.canvasRef.current.getContext("2d");

    ctx.clearRect(0, 0, width, height);

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
    const { obsOrVarContinuousFieldDisplayName, bins } = this.props;
    if (
      prevProps.obsOrVarContinuousFieldDisplayName !==
        obsOrVarContinuousFieldDisplayName ||
      prevProps.bins !== bins
    )
      this.drawHistogram();
  };

  render() {
    const {
      domainLabel,
      obsOrVarContinuousFieldDisplayName,
      width,
      height,
    } = this.props;

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
          style={{
            marginRight: 5,
            width,
            height,
            borderBottom: "solid rgb(230, 230, 230) 0.25px",
          }}
          className="mini-histo"
          width={width}
          height={height}
          ref={this.canvasRef}
        />
        <div key="text" style={{ fontSize: "14px" }}>
          <p style={{ margin: "0" }}>
            This histograms shows the distribution of{" "}
            <strong>{obsOrVarContinuousFieldDisplayName}</strong> within{" "}
            <strong>{domainLabel}</strong>.
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
