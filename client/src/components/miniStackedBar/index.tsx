import React from "react";

export default class MiniStackedBar extends React.PureComponent {
  canvasRef: any;

  constructor(props: any) {
    super(props);
    this.canvasRef = React.createRef();
  }

  drawStacks = () => {
    const {
      // @ts-expect-error ts-migrate(2339) FIXME: Property 'domainValues' does not exist on type 'Re... Remove this comment to see the full error message
      domainValues,
      // @ts-expect-error ts-migrate(2339) FIXME: Property 'scale' does not exist on type 'Readonly<... Remove this comment to see the full error message
      scale,
      // @ts-expect-error ts-migrate(2339) FIXME: Property 'domain' does not exist on type 'Readonly... Remove this comment to see the full error message
      domain,
      // @ts-expect-error ts-migrate(2339) FIXME: Property 'colorTable' does not exist on type 'Read... Remove this comment to see the full error message
      colorTable,
      // @ts-expect-error ts-migrate(2339) FIXME: Property 'occupancy' does not exist on type 'Reado... Remove this comment to see the full error message
      occupancy,
      // @ts-expect-error ts-migrate(2339) FIXME: Property 'width' does not exist on type 'Readonly<... Remove this comment to see the full error message
      width,
      // @ts-expect-error ts-migrate(2339) FIXME: Property 'height' does not exist on type 'Readonly... Remove this comment to see the full error message
      height,
    } = this.props;

    if (!colorTable || !domainValues) return;

    const { scale: colorScale } = colorTable;

    const ctx = this.canvasRef?.current.getContext("2d");

    ctx.clearRect(0, 0, width, height);
    let currentOffset = 0;

    let occupancyValue;
    let scaledValue;
    let value;

    for (let i = 0, { length } = domainValues; i < length; i += 1) {
      value = domainValues[i];
      occupancyValue = occupancy.get(value);
      scaledValue = scale(occupancyValue);
      ctx.fillStyle = occupancyValue
        ? colorScale(domain.indexOf(value))
        : "rgb(255,255,255)";
      ctx.fillRect(currentOffset, 0, occupancyValue ? scaledValue : 0, height);
      currentOffset += occupancyValue ? scaledValue : 0;
    }
  };

  componentDidUpdate = (prevProps: any) => {
    // @ts-expect-error ts-migrate(2339) FIXME: Property 'occupancy' does not exist on type 'Reado... Remove this comment to see the full error message
    const { occupancy } = this.props;
    if (occupancy !== prevProps.occupancy) this.drawStacks();
  };

  componentDidMount = () => {
    this.drawStacks();
  };

  render() {
    // @ts-expect-error ts-migrate(2339) FIXME: Property 'width' does not exist on type 'Readonly<... Remove this comment to see the full error message
    const { width, height } = this.props;
    // @ts-expect-error ts-migrate(2339) FIXME: Property 'canvas' does not exist on type 'MiniStac... Remove this comment to see the full error message
    const { canvas } = this;
    if (canvas) canvas.getContext("2d").clearRect(0, 0, width, height);

    return (
      <canvas
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
