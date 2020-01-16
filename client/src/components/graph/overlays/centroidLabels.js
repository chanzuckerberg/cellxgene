/* eslint-disable jsx-a11y/mouse-events-have-key-events */
import React, { PureComponent } from "react";
import { connect } from "react-redux";

export default
@connect(state => ({
  colorAccessor: state.colors.colorAccessor,
  dilatedValue: state.pointDilation.categoryField,
  labels: state.centroidLabels.labels
}))
class CentroidLabels extends PureComponent {
  render() {
    const {
      labels,
      inverseTransform,
      dilatedValue,
      dispatch,
      colorAccessor
    } = this.props;

    const labelSVGS = [];
    let fontSize = "15px";
    let fontWeight = null;
    labels.forEach((value, key) => {
      fontSize = "15px";
      fontWeight = null;
      if (key === dilatedValue) {
        fontSize = "18px";
        fontWeight = "800";
      }
      labelSVGS.push(
        <g
          // eslint-disable-next-line react/no-array-index-key
          key={key}
          className="centroid-label"
          transform={`translate(${value[0]}, ${value[1]})`}
        >
          <text
            transform={inverseTransform}
            textAnchor="middle"
            data-label={key}
            style={{
              fontFamily: "Roboto Condensed",
              fontSize,
              fontWeight,
              fill: "black",
              userSelect: "none"
            }}
            onMouseEnter={e =>
              dispatch({
                type: "category value mouse hover start",
                metadataField: colorAccessor,
                categoryField: e.target.getAttribute("data-label")
              })
            }
            onMouseOut={e =>
              dispatch({
                type: "category value mouse hover end",
                metadataField: colorAccessor,
                categoryField: e.target.getAttribute("data-label")
              })
            }
            pointerEvents="visiblePainted"
          >
            {key.length > 20 ? `${key.substr(0, 20)}...` : key}
          </text>
        </g>
      );
    });

    return <>{labelSVGS}</>;
  }
}
