import React from "react";
import {
  Popover,
  PopoverInteractionKind,
  Position,
  Classes,
} from "@blueprintjs/core";

export default class MiniHistogram extends React.PureComponent {
  // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
  canvasRef: any;

  // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any -- - FIXME: disabled temporarily on migrate to TS.
  constructor(props: any) {
    super(props);
    this.canvasRef = React.createRef();
  }

  // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types --- FIXME: disabled temporarily on migrate to TS.
  drawHistogram = () => {
    // @ts-expect-error ts-migrate(2339) FIXME: Property 'xScale' does not exist on type 'Readonly... Remove this comment to see the full error message
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

  // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types --- FIXME: disabled temporarily on migrate to TS.
  componentDidMount = () => {
    this.drawHistogram();
  };

  // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any -- - FIXME: disabled temporarily on migrate to TS.
  componentDidUpdate = (prevProps: any) => {
    // @ts-expect-error ts-migrate(2339) FIXME: Property 'obsOrVarContinuousFieldDisplayName' does... Remove this comment to see the full error message
    const { obsOrVarContinuousFieldDisplayName, bins } = this.props;
    if (
      prevProps.obsOrVarContinuousFieldDisplayName !==
        obsOrVarContinuousFieldDisplayName ||
      prevProps.bins !== bins
    )
      this.drawHistogram();
  };

  // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types --- FIXME: disabled temporarily on migrate to TS.
  render() {
    const {
      // @ts-expect-error ts-migrate(2339) FIXME: Property 'domainLabel' does not exist on type 'Rea... Remove this comment to see the full error message
      domainLabel,
      // @ts-expect-error ts-migrate(2339) FIXME: Property 'obsOrVarContinuousFieldDisplayName' does... Remove this comment to see the full error message
      obsOrVarContinuousFieldDisplayName,
      // @ts-expect-error ts-migrate(2339) FIXME: Property 'width' does not exist on type 'Readonly<... Remove this comment to see the full error message
      width,
      // @ts-expect-error ts-migrate(2339) FIXME: Property 'height' does not exist on type 'Readonly... Remove this comment to see the full error message
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
