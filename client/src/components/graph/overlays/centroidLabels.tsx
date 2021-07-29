import React, { PureComponent } from "react";
import { connect, shallowEqual } from "react-redux";
import Async from "react-async";

import { categoryLabelDisplayStringLongLength } from "../../../globals";
import calcCentroid from "../../../util/centroid";
import { createColorQuery } from "../../../util/stateManager/colorHelpers";

// @ts-expect-error ts-migrate(1238) FIXME: Unable to resolve signature of class decorator whe... Remove this comment to see the full error message
@connect((state) => ({
  // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
  annoMatrix: (state as any).annoMatrix,
  // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
  colors: (state as any).colors,
  // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
  layoutChoice: (state as any).layoutChoice,
  // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
  dilatedValue: (state as any).pointDilation.categoryField,
  // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
  categoricalSelection: (state as any).categoricalSelection,
  // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
  showLabels: (state as any).centroidLabels?.showLabels,
  // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
  genesets: (state as any).genesets.genesets,
}))
export default class CentroidLabels extends PureComponent {
  // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any -- - FIXME: disabled temporarily on migrate to TS.
  static watchAsync(props: any, prevProps: any) {
    return !shallowEqual(props.watchProps, prevProps.watchProps);
  }

  // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any -- - FIXME: disabled temporarily on migrate to TS.
  fetchAsyncProps = async (props: any) => {
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
    // @ts-expect-error ts-migrate(2339) FIXME: Property 'overlaySetShowing' does not exist on typ... Remove this comment to see the full error message
    const { overlaySetShowing } = this.props;
    overlaySetShowing("centroidLabels", showLabels && labels.size > 0);
    return {
      labels,
      colorAccessor,
      category: categoricalSelection[colorAccessor],
    };
  };

  // @ts-expect-error ts-migrate(6133) FIXME: 'e' is declared but its value is never read.
  // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any -- - FIXME: disabled temporarily on migrate to TS.
  handleMouseEnter = (e: any, colorAccessor: any, label: any) => {
    // @ts-expect-error ts-migrate(2339) FIXME: Property 'dispatch' does not exist on type 'Readon... Remove this comment to see the full error message
    const { dispatch } = this.props;
    dispatch({
      type: "category value mouse hover start",
      metadataField: colorAccessor,
      categoryField: label,
    });
  };

  // @ts-expect-error ts-migrate(6133) FIXME: 'e' is declared but its value is never read.
  // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any -- - FIXME: disabled temporarily on migrate to TS.
  handleMouseOut = (e: any, colorAccessor: any, label: any) => {
    // @ts-expect-error ts-migrate(2339) FIXME: Property 'dispatch' does not exist on type 'Readon... Remove this comment to see the full error message
    const { dispatch } = this.props;
    dispatch({
      type: "category value mouse hover end",
      metadataField: colorAccessor,
      categoryField: label,
    });
  };

  // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types --- FIXME: disabled temporarily on migrate to TS.
  colorByQuery() {
    // @ts-expect-error ts-migrate(2339) FIXME: Property 'annoMatrix' does not exist on type 'Read... Remove this comment to see the full error message
    const { annoMatrix, colors, genesets } = this.props;
    const { schema } = annoMatrix;
    const { colorMode, colorAccessor } = colors;
    return createColorQuery(colorMode, colorAccessor, schema, genesets);
  }

  // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types --- FIXME: disabled temporarily on migrate to TS.
  async fetchData() {
    // @ts-expect-error ts-migrate(2339) FIXME: Property 'annoMatrix' does not exist on type 'Read... Remove this comment to see the full error message
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

  // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types --- FIXME: disabled temporarily on migrate to TS.
  render() {
    const {
      // @ts-expect-error ts-migrate(2339) FIXME: Property 'inverseTransform' does not exist on type... Remove this comment to see the full error message
      inverseTransform,
      // @ts-expect-error ts-migrate(2339) FIXME: Property 'dilatedValue' does not exist on type 'Re... Remove this comment to see the full error message
      dilatedValue,
      // @ts-expect-error ts-migrate(2339) FIXME: Property 'categoricalSelection' does not exist on ... Remove this comment to see the full error message
      categoricalSelection,
      // @ts-expect-error ts-migrate(2339) FIXME: Property 'showLabels' does not exist on type 'Read... Remove this comment to see the full error message
      showLabels,
      // @ts-expect-error ts-migrate(2339) FIXME: Property 'colors' does not exist on type 'Readonly... Remove this comment to see the full error message
      colors,
      // @ts-expect-error ts-migrate(2339) FIXME: Property 'annoMatrix' does not exist on type 'Read... Remove this comment to see the full error message
      annoMatrix,
      // @ts-expect-error ts-migrate(2339) FIXME: Property 'layoutChoice' does not exist on type 'Re... Remove this comment to see the full error message
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
            // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
            const labelSVGS: any = [];
            const deselectOpacity = 0.375;
            // @ts-expect-error ts-migrate(2339) FIXME: Property 'category' does not exist on type 'unknow... Remove this comment to see the full error message
            const { category, colorAccessor, labels } = asyncProps;
            // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
            labels.forEach((coords: any, label: any) => {
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
  onMouseOut, // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
}: any) => {
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
          // @ts-expect-error ts-migrate(2322) FIXME: Type 'string | null' is not assignable to type 'Fo... Remove this comment to see the full error message
          fontWeight,
          fill: "black",
          userSelect: "none",
          // @ts-expect-error ts-migrate(2322) FIXME: Type '{ opacity: any; }' is not assignable to type... Remove this comment to see the full error message
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
