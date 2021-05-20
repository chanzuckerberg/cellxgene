import React, { PureComponent } from "react";
import { connect, shallowEqual } from "react-redux";
import Async from "react-async";

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
  genesets: state.genesets.genesets,
  differential: state.differential,
}))
class CentroidLabels extends PureComponent {
  static watchAsync(props, prevProps) {
    return !shallowEqual(props.watchProps, prevProps.watchProps);
  }

  fetchAsyncProps = async (props) => {
    const {
      annoMatrix,
      colors,
      layoutChoice,
      categoricalSelection,
      showLabels,
    } = props.watchProps;
    const { schema } = annoMatrix;
    const { colorAccessor } = colors;

    const [layoutDf, colorDf] = await this.fetchData();
    let labels;
    if (colorDf) {
      labels = calcCentroid(
        schema,
        colorAccessor,
        colorDf,
        layoutChoice,
        layoutDf
      );
    } else {
      labels = new Map();
    }

    const { overlaySetShowing } = this.props;
    overlaySetShowing("centroidLabels", showLabels && labels.size > 0);

    return {
      labels,
      colorAccessor,
      category: categoricalSelection[colorAccessor],
    };
  };

  handleMouseEnter = (e, colorAccessor, label) => {
    const { dispatch } = this.props;
    dispatch({
      type: "category value mouse hover start",
      metadataField: colorAccessor,
      categoryField: label,
    });
  };

  handleMouseOut = (e, colorAccessor, label) => {
    const { dispatch } = this.props;
    dispatch({
      type: "category value mouse hover end",
      metadataField: colorAccessor,
      categoryField: label,
    });
  };

  colorByQuery() {
    const { annoMatrix, colors, genesets } = this.props;
    const { schema } = annoMatrix;
    const { colorMode, colorAccessor } = colors;
    return createColorQuery(colorMode, colorAccessor, schema, genesets);
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

  render() {
    const {
      inverseTransform,
      dilatedValue,
      categoricalSelection,
      showLabels,
      colors,
      annoMatrix,
      layoutChoice,
    } = this.props;

    return (
      <Async
        watchFn={CentroidLabels.watchAsync}
        promiseFn={this.fetchAsyncProps}
        watchProps={{
          annoMatrix,
          colors,
          layoutChoice,
          categoricalSelection,
          dilatedValue,
          showLabels,
        }}
      >
        <Async.Fulfilled>
          {(asyncProps) => {
            if (!showLabels) return null;

            const labelSVGS = [];
            const deselectOpacity = 0.375;
            const { category, colorAccessor, labels } = asyncProps;

            labels.forEach((coords, label) => {
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
                // eslint-disable-next-line jsx-a11y/mouse-events-have-key-events -- the mouse actions for centroid labels do not have a screen reader alternative
                <Label
                  key={label} // eslint-disable-line react/no-array-index-key --- label is not an index, eslint is confused
                  label={label}
                  dilatedValue={dilatedValue}
                  coords={coords}
                  inverseTransform={inverseTransform}
                  opactity={selected ? 1 : deselectOpacity}
                  colorAccessor={colorAccessor}
                  displayLabel={displayLabel}
                  onMouseEnter={this.handleMouseEnter}
                  onMouseOut={this.handleMouseOut}
                />
              );
            });

            return <>{labelSVGS}</>;
          }}
        </Async.Fulfilled>
      </Async>
    );
  }
}

const Label = ({
  label,
  dilatedValue,
  coords,
  inverseTransform,
  opacity,
  colorAccessor,
  displayLabel,
  onMouseEnter,
  onMouseOut,
}) => {
  /*
  Render a label at a given coordinate.
  */
  let fontSize = "15px";
  let fontWeight = null;
  if (label === dilatedValue) {
    fontSize = "18px";
    fontWeight = "800";
  }

  return (
    <g
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
        style={{
          fontSize,
          fontWeight,
          fill: "black",
          userSelect: "none",
          opacity: { opacity },
        }}
        onMouseEnter={(e) => onMouseEnter(e, colorAccessor, label)}
        onMouseOut={(e) => onMouseOut(e, colorAccessor, label)}
        pointerEvents="visiblePainted"
      >
        {displayLabel}
      </text>
    </g>
  );
};
