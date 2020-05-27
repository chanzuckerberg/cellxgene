import React, { PureComponent } from "react";
import { connect } from "react-redux";

import { categoryLabelDisplayStringLongLength } from "../../../globals";

export default
@connect((state) => ({
  colorAccessor: state.colors.colorAccessor,
  dilatedValue: state.pointDilation.categoryField,
  labels: state.centroidLabels.labels,
  categoricalSelection: state.categoricalSelection,
}))
class CentroidLabels extends PureComponent {
  // Check to see if centroids have either just been displayed or removed from the overlay

  componentDidUpdate(prevProps) {
    const { labels, overlayToggled } = this.props;
    const prevSize = prevProps.labels.size;
    const { size } = labels;

    const displayChangeOff = prevSize > 0 && size === undefined;
    const displayChangeOn = prevSize === undefined && size > 0;

    if (displayChangeOn || displayChangeOff) {
      // Notify overlay layer of display change
      overlayToggled("centroidLabels", displayChangeOn);
    }
  }

  render() {
    const {
      labels,
      inverseTransform,
      dilatedValue,
      dispatch,
      colorAccessor,
      categoricalSelection,
    } = this.props;

    if (!colorAccessor || labels.size === undefined || labels.size === 0)
      return null;

    const {
      categoryValueIndices,
      categoryValueSelected,
    } = categoricalSelection?.[colorAccessor];

    const labelSVGS = [];
    let fontSize = "15px";
    let fontWeight = null;
    const deselectOpacity = 0.375;
    labels.forEach((coords, label) => {
      fontSize = "15px";
      fontWeight = null;
      if (label === dilatedValue) {
        fontSize = "18px";
        fontWeight = "800";
      }

      const selected = categoryValueSelected[categoryValueIndices.get(label)];

      // Mirror LSB middle truncation
      let displayLabel = label;
      if (displayLabel.length > categoryLabelDisplayStringLongLength) {
        displayLabel = `${label.slice(
          0,
          categoryLabelDisplayStringLongLength / 2
        )}…${label.slice(-categoryLabelDisplayStringLongLength / 2)}`;
      }

      labelSVGS.push(
        <g
          // label is unique so disabling eslint rule
          // eslint-disable-next-line react/no-array-index-key
          key={label}
          className="centroid-label"
          transform={`translate(${coords[0]}, ${coords[1]})`}
          data-testclass="centroid-label"
          data-testid={`${label}-centroid-label`}
        >
          {/* The mouse actions for centroid labels do not have a screen reader alternative */}
          {/* eslint-disable-next-line jsx-a11y/mouse-events-have-key-events */}
          <text
            transform={inverseTransform}
            textAnchor="middle"
            data-label={label}
            style={{
              fontSize,
              fontWeight,
              fill: "black",
              userSelect: "none",
              opacity: selected ? 1 : deselectOpacity,
            }}
            onMouseEnter={(e) =>
              dispatch({
                type: "category value mouse hover start",
                metadataField: colorAccessor,
                categoryField: e.target.getAttribute("data-label"),
              })
            }
            onMouseOut={(e) =>
              dispatch({
                type: "category value mouse hover end",
                metadataField: colorAccessor,
                categoryField: e.target.getAttribute("data-label"),
              })
            }
            pointerEvents="visiblePainted"
          >
            {displayLabel}
          </text>
        </g>
      );
    });

    return <>{labelSVGS}</>;
  }
}
