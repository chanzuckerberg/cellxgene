/* eslint-disable jsx-a11y/mouse-events-have-key-events */
import React, { PureComponent } from "react";

export default class CentroidLabels extends PureComponent {
  render() {
    const fontSize = "20px";

    const { labels, inverseScale, mouseEnter, mouseExit } = this.props;

    const iter = labels.entries();
    let pair = iter.next().value;
    let value;
    let key;

    const labelSVGS = [];
    while (pair) {
      key = pair[0];
      value = pair[1];
      labelSVGS.push(
        <g
          key={key}
          className="centroid-label"
          transform={`translate(${value[0]}, ${value[1]})`}
        >
          <text
            transform={inverseScale}
            textAnchor="middle"
            data-label={key}
            id={`svg${key.replace(/[^\w]/gi, "")}-label`}
            style={{
              fontFamily: "Roboto Condensed",
              fontSize,
              fill: "black"
            }}
            onMouseEnter={mouseEnter}
            onMouseOut={mouseExit}
            pointerEvents="visiblePainted"
          >
            {key.length > 20 ? `${key.substr(0, 20)}...` : key}
          </text>
        </g>
      );
      pair = iter.next().value;
    }

    return <>{labelSVGS}</>;
  }
}
