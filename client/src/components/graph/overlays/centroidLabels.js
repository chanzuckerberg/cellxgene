import React, { PureComponent } from "react";
import { connect } from "react-redux";

import { categoryLabelDisplayStringLongLength } from "../../../globals";
import calcCentroid from "../../../util/centroid";
import { createColorQuery } from "../../../util/stateManager/colorHelpers";

export default
@connect((state) => ({
  annoMatrix: state.annoMatrix,
  colors: state.colors,
  layoutChoice: state.layoutChoice,
  dilatedValue: state.pointDilation.categoryField,
  categoricalSelection: state.categoricalSelection,
  showLabels: state.centroidLabels?.showLabels,
}))
class CentroidLabels extends PureComponent {
  constructor(props) {
    super(props);
    this.state = {
      status: "pending",
      labels: new Map(),
      colorAccessor: null,
    };
  }

  colorByQuery() {
    const { annoMatrix, colors } = this.props;
    const { schema } = annoMatrix;
    const { colorMode, colorAccessor } = colors;
    return createColorQuery(colorMode, colorAccessor, schema);
  }

  async fetchData() {
    const { annoMatrix, layoutChoice } = this.props;
    // fetch all data we need: layout, category
    const promises = [];
    // layout
    promises.push(annoMatrix.fetch("emb", layoutChoice.current));
    // category to label - we ONLY label on obs, never on X, etc.
    const query = this.colorByQuery();
    if (query && query[0] === "obs") {
      promises.push(annoMatrix.fetch(...query));
    } else {
      promises.push(Promise.resolve(null));
    }

    return Promise.all(promises);
  }

  // Check to see if centroids have either just been displayed or removed from the overlay
  async updateState(prevProps) {
    const { annoMatrix, colors, layoutChoice } = this.props;
    if (!annoMatrix || !layoutChoice || !colors) return;

    if (
      annoMatrix !== prevProps?.annoMatrix ||
      colors !== prevProps?.colors ||
      layoutChoice !== prevProps?.layoutChoice
    ) {
      const { schema } = annoMatrix;
      const { colorAccessor } = colors;
      this.setState({ status: "pending" });
      try {
        const [layoutDf, colorDf] = await this.fetchData();
        let labels = [];
        if (colorDf) {
          labels = calcCentroid(
            schema,
            colorAccessor,
            colorDf,
            layoutChoice,
            layoutDf
          );
        }
        this.setState({ status: "success", labels, colorAccessor });
      } catch (error) {
        this.setState({ status: "error", error });
        throw error;
      }
    }
    return;
  }

  componentDidMount() {
    this.updateState(null);
  }

  componentDidUpdate(prevProps, prevState) {
    this.updateState(prevProps);

    const { showLabels, overlaySetShowing } = this.props;
    const { labels } = this.state;
    overlaySetShowing("centroidLabels", showLabels && labels.size > 0);
  }

  render() {
    const {
      inverseTransform,
      dilatedValue,
      dispatch,
      colors,
      categoricalSelection,
      showLabels,
    } = this.props;

    const { status, labels, colorAccessor } = this.state;

    if (
      status !== "success" ||
      !showLabels ||
      !colorAccessor ||
      labels.size === 0
    )
      return null;

    const category = categoricalSelection[colorAccessor];

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

      const selected = category.get(label) ?? true;

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
          // eslint-disable-next-line react/no-array-index-key --- label is unique and consistent
          key={label}
          className="centroid-label"
          transform={`translate(${coords[0]}, ${coords[1]})`}
          data-testclass="centroid-label"
          data-testid={`${label}-centroid-label`}
        >
          {/* eslint-disable-next-line jsx-a11y/mouse-events-have-key-events --- the mouse actions for centroid labels do not have a screen reader alternative*/}
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
