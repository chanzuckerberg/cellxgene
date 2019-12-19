/* eslint-disable jsx-a11y/mouse-events-have-key-events */
import React, { PureComponent } from "react";

export default class CentroidLabels extends PureComponent {
  render() {
    const {
      labels,
      inverseTransform,
      mouseEnter,
      mouseExit,
      dilatedValue
    } = this.props;

    const iter = labels.entries();
    let pair = iter.next().value;
    let value;
    let key;

    const labelSVGS = [];
    let fontSize = "15px";
    let fontWeight = null;
    while (pair) {
      key = pair[0];
      value = pair[1];
      fontSize = "15px";
      fontWeight = null;
      if (key === dilatedValue) {
        fontSize = "18px";
        fontWeight = "800";
      }
      labelSVGS.push(
        <g
          key={key}
          className="centroid-label"
          transform={`translate(${value[0]}, ${value[1]})`}
        >
          <text
            transform={inverseTransform}
            textAnchor="middle"
            data-label={key}
            id={`svg${key.replace(/[^\w]/gi, "")}-label`}
            style={{
              fontFamily: "Roboto Condensed",
              fontSize,
              fontWeight,
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
