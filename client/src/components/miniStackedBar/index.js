import React from "react";

export default class MiniStackedBar extends React.PureComponent {
  constructor(props) {
    super(props);
    this.canvasRef = React.createRef();
  }

  drawStacks = () => {
    const {
      domainValues,
      scale,
      domain,
      colorTable,
      occupancy,
      width,
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

  componentDidUpdate = (prevProps) => {
    const { occupancy } = this.props;
    if (occupancy !== prevProps.occupancy) this.drawStacks();
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
