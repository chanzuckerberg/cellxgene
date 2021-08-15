import React from "react";
import { connect, shallowEqual } from "react-redux";
import * as d3 from "d3";
import Async from "react-async";
import memoize from "memoize-one";
import * as globals from "../../globals";
import actions from "../../actions";
import { makeContinuousDimensionName } from "../../util/nameCreators";
import HistogramHeader from "./header";
import Histogram from "./histogram";
import HistogramFooter from "./footer";
import StillLoading from "./loading";
import ErrorLoading from "./error";
import { Dataframe } from "../../util/dataframe";

const MARGIN = {
  LEFT: 10, // Space for 0 tick label on X axis
  RIGHT: 54, // space for Y axis & labels
  BOTTOM: 25, // space for X axis & labels
  TOP: 3,
};
const WIDTH = 340 - MARGIN.LEFT - MARGIN.RIGHT;
const HEIGHT = 135 - MARGIN.TOP - MARGIN.BOTTOM;
const MARGIN_MINI = {
  LEFT: 0, // Space for 0 tick label on X axis
  RIGHT: 0, // space for Y axis & labels
  BOTTOM: 0, // space for X axis & labels
  TOP: 0,
};
const WIDTH_MINI = 120 - MARGIN_MINI.LEFT - MARGIN_MINI.RIGHT;
const HEIGHT_MINI = 15 - MARGIN_MINI.TOP - MARGIN_MINI.BOTTOM;

// @ts-expect-error ts-migrate(1238) FIXME: Unable to resolve signature of class decorator whe... Remove this comment to see the full error message
@connect((state, ownProps) => {
  // @ts-expect-error ts-migrate(2339) FIXME: Property 'isObs' does not exist on type '{}'.
  const { isObs, isUserDefined, isGeneSetSummary, field } = ownProps;
  const myName = makeContinuousDimensionName(
    { isObs, isUserDefined, isGeneSetSummary },
    field
  );
  return {
    // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
    annoMatrix: (state as any).annoMatrix,
    isScatterplotXXaccessor:
      // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
      (state as any).controls.scatterplotXXaccessor === field,
    isScatterplotYYaccessor:
      // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
      (state as any).controls.scatterplotYYaccessor === field,
    // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
    continuousSelectionRange: (state as any).continuousSelection[myName],
    // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
    isColorAccessor: (state as any).colors.colorAccessor === field,
  };
})
class HistogramBrush extends React.PureComponent {
  // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any -- - FIXME: disabled temporarily on migrate to TS.
  static watchAsync(props: any, prevProps: any) {
    return !shallowEqual(props.watchProps, prevProps.watchProps);
  }

  /* memoized closure to prevent HistogramHeader unecessary repaint */
  // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
  handleColorAction = memoize((dispatch) => (field: any, isObs: any) => {
    if (isObs) {
      dispatch({
        type: "color by continuous metadata",
        colorAccessor: field,
      });
    } else {
      dispatch(actions.requestSingleGeneExpressionCountsForColoringPOST(field));
    }
  });

  // @ts-expect-error ts-migrate(6133) FIXME: 'selection' is declared but its value is never rea... Remove this comment to see the full error message
  // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any -- - FIXME: disabled temporarily on migrate to TS.
  onBrush = (selection: any, x: any, eventType: any) => {
    const type = `continuous metadata histogram ${eventType}`;
    return () => {
      const {
        // @ts-expect-error ts-migrate(2339) FIXME: Property 'dispatch' does not exist on type 'Readon... Remove this comment to see the full error message
        dispatch,
        // @ts-expect-error ts-migrate(2339) FIXME: Property 'field' does not exist on type 'Readonly<... Remove this comment to see the full error message
        field,
        // @ts-expect-error ts-migrate(2339) FIXME: Property 'isObs' does not exist on type 'Readonly<... Remove this comment to see the full error message
        isObs,
        // @ts-expect-error ts-migrate(2339) FIXME: Property 'isUserDefined' does not exist on type 'R... Remove this comment to see the full error message
        isUserDefined,
        // @ts-expect-error ts-migrate(2339) FIXME: Property 'isGeneSetSummary' does not exist on type... Remove this comment to see the full error message
        isGeneSetSummary,
      } = this.props;

      // ignore programmatically generated events
      // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
      if (!(d3 as any).event.sourceEvent) return;
      // ignore cascading events, which are programmatically generated
      // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
      if ((d3 as any).event.sourceEvent.sourceEvent) return;

      const query = this.createQuery();
      // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
      const range = (d3 as any).event.selection
        ? // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
          [x((d3 as any).event.selection[0]), x((d3 as any).event.selection[1])]
        : null;
      const otherProps = {
        selection: field,
        continuousNamespace: {
          isObs,
          isUserDefined,
          isGeneSetSummary,
        },
      };
      dispatch(
        actions.selectContinuousMetadataAction(type, query, range, otherProps)
      );
    };
  };

  // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any -- - FIXME: disabled temporarily on migrate to TS.
  onBrushEnd =
    (
      _selection: any, // eslint-disable-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any -- - FIXME: disabled temporarily on migrate to TS.
      x: any // eslint-disable-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any -- - FIXME: disabled temporarily on migrate to TS.
    ) =>
    // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types --- FIXME: disabled temporarily on migrate to TS.
    () => {
      const {
        // @ts-expect-error ts-migrate(2339) FIXME: Property 'dispatch' does not exist on type 'Readon... Remove this comment to see the full error message
        dispatch,
        // @ts-expect-error ts-migrate(2339) FIXME: Property 'field' does not exist on type 'Readonly<... Remove this comment to see the full error message
        field,
        // @ts-expect-error ts-migrate(2339) FIXME: Property 'isObs' does not exist on type 'Readonly<... Remove this comment to see the full error message
        isObs,
        // @ts-expect-error ts-migrate(2339) FIXME: Property 'isUserDefined' does not exist on type 'R... Remove this comment to see the full error message
        isUserDefined,
        // @ts-expect-error ts-migrate(2339) FIXME: Property 'isGeneSetSummary' does not exist on type... Remove this comment to see the full error message
        isGeneSetSummary,
      } = this.props;
      const minAllowedBrushSize = 10;
      const smallAmountToAvoidInfiniteLoop = 0.1;

      // ignore programmatically generated events
      // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
      if (!(d3 as any).event.sourceEvent) return;
      // ignore cascading events, which are programmatically generated
      // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
      if ((d3 as any).event.sourceEvent.sourceEvent) return;

      let type;
      let range = null;
      // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
      if ((d3 as any).event.selection) {
        type = "continuous metadata histogram end";
        if (
          // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
          (d3 as any).event.selection[1] - (d3 as any).event.selection[0] >
          minAllowedBrushSize
        ) {
          range = [
            // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
            x((d3 as any).event.selection[0]),
            // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
            x((d3 as any).event.selection[1]),
          ];
        } else {
          /* the user selected range is too small and will be hidden #587, so take control of it procedurally */
          /* https://stackoverflow.com/questions/12354729/d3-js-limit-size-of-brush */
          const procedurallyResizedBrushWidth =
            // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
            (d3 as any).event.selection[0] +
            minAllowedBrushSize +
            smallAmountToAvoidInfiniteLoop; //
          range = [
            // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
            x((d3 as any).event.selection[0]),
            x(procedurallyResizedBrushWidth),
          ];
        }
      } else {
        type = "continuous metadata histogram cancel";
      }

      const query = this.createQuery();
      const otherProps = {
        selection: field,
        continuousNamespace: {
          isObs,
          isUserDefined,
          isGeneSetSummary,
        },
      };
      dispatch(
        actions.selectContinuousMetadataAction(type, query, range, otherProps)
      );
    };

  // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types --- FIXME: disabled temporarily on migrate to TS.
  handleSetGeneAsScatterplotX = () => {
    // @ts-expect-error ts-migrate(2339) FIXME: Property 'dispatch' does not exist on type 'Readon... Remove this comment to see the full error message
    const { dispatch, field } = this.props;
    dispatch({
      type: "set scatterplot x",
      data: field,
    });
  };

  // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types --- FIXME: disabled temporarily on migrate to TS.
  handleSetGeneAsScatterplotY = () => {
    // @ts-expect-error ts-migrate(2339) FIXME: Property 'dispatch' does not exist on type 'Readon... Remove this comment to see the full error message
    const { dispatch, field } = this.props;
    dispatch({
      type: "set scatterplot y",
      data: field,
    });
  };

  // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types --- FIXME: disabled temporarily on migrate to TS.
  removeHistogram = () => {
    const {
      // @ts-expect-error ts-migrate(2339) FIXME: Property 'dispatch' does not exist on type 'Readon... Remove this comment to see the full error message
      dispatch,
      // @ts-expect-error ts-migrate(2339) FIXME: Property 'field' does not exist on type 'Readonly<... Remove this comment to see the full error message
      field,
      // @ts-expect-error ts-migrate(2339) FIXME: Property 'isColorAccessor' does not exist on type ... Remove this comment to see the full error message
      isColorAccessor,
      // @ts-expect-error ts-migrate(2339) FIXME: Property 'isScatterplotXXaccessor' does not exist ... Remove this comment to see the full error message
      isScatterplotXXaccessor,
      // @ts-expect-error ts-migrate(2339) FIXME: Property 'isScatterplotYYaccessor' does not exist ... Remove this comment to see the full error message
      isScatterplotYYaccessor,
    } = this.props;
    dispatch({
      type: "clear user defined gene",
      data: field,
    });
    if (isColorAccessor) {
      dispatch({
        type: "reset colorscale",
      });
    }
    if (isScatterplotXXaccessor) {
      dispatch({
        type: "set scatterplot x",
        data: null,
      });
    }
    if (isScatterplotYYaccessor) {
      dispatch({
        type: "set scatterplot y",
        data: null,
      });
    }
  };

  // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types --- FIXME: disabled temporarily on migrate to TS.
  fetchAsyncProps = async () => {
    // @ts-expect-error ts-migrate(2339) FIXME: Property 'annoMatrix' does not exist on type 'Read... Remove this comment to see the full error message
    const { annoMatrix, width } = this.props;

    const { isClipped } = annoMatrix;

    const query = this.createQuery();
    // @ts-expect-error ts-migrate(2488) FIXME: Type 'any[] | null' must have a '[Symbol.iterator]... Remove this comment to see the full error message
    const df: Dataframe = await annoMatrix.fetch(...query);
    const column = df.icol(0);

    // if we are clipped, fetch both our value and our unclipped value,
    // as we need the absolute min/max range, not just the clipped min/max.
    const summary = column.summarizeContinuous();
    const range = [summary.min, summary.max];

    let unclippedRange = [...range];
    if (isClipped) {
      const parent: Dataframe = await annoMatrix.viewOf.fetch(...query);
      const { min, max } = parent.icol(0).summarizeContinuous();
      unclippedRange = [min, max];
    }

    const unclippedRangeColor = [
      !annoMatrix.isClipped || annoMatrix.clipRange[0] === 0
        ? "#bbb"
        : globals.blue,
      !annoMatrix.isClipped || annoMatrix.clipRange[1] === 1
        ? "#bbb"
        : globals.blue,
    ];

    const histogram = this.calcHistogramCache(
      column,
      MARGIN,
      width || WIDTH,
      HEIGHT
    );
    const miniHistogram = this.calcHistogramCache(
      column,
      MARGIN_MINI,
      width || WIDTH_MINI,
      HEIGHT_MINI
    );

    const isSingleValue = summary.min === summary.max;
    const nonFiniteExtent =
      summary.min === undefined ||
      summary.max === undefined ||
      Number.isNaN(summary.min) ||
      Number.isNaN(summary.max);

    const OK2Render = !summary.categorical && !nonFiniteExtent;

    return {
      histogram,
      miniHistogram,
      range,
      unclippedRange,
      unclippedRangeColor,
      isSingleValue,
      OK2Render,
    };
  };

  // eslint-disable-next-line class-methods-use-this, @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any -- instance method allows for memoization per annotation
  calcHistogramCache(col: any, newMargin: any, newWidth: any, newHeight: any) {
    /*
     recalculate expensive stuff, notably bins, summaries, etc.
    */
    const histogramCache = {}; /* maybe change this so that it computes ... */
    const summary =
      col.summarizeContinuous(); /* this is memoized, so it's free the second time you call it */
    const { min: domainMin, max: domainMax } = summary;
    const numBins = 40;
    const { TOP: topMargin, LEFT: leftMargin } = newMargin;
    // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
    (histogramCache as any).domain = [domainMin, domainMax];
    // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
    /* doesn't change with mini */ (histogramCache as any).x = d3
      .scaleLinear()
      .domain([domainMin, domainMax])
      .range([leftMargin, leftMargin + newWidth]);

    // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
    (histogramCache as any).bins = col.histogramContinuous(numBins, [
      domainMin,
      domainMax,
    ]);
    // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
    /* memoized */ (histogramCache as any).binWidth =
      (domainMax - domainMin) / numBins;

    // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
    (histogramCache as any).binStart = (i: any) =>
      // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
      domainMin + i * (histogramCache as any).binWidth;
    // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
    (histogramCache as any).binEnd = (i: any) =>
      // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
      domainMin + (i + 1) * (histogramCache as any).binWidth;

    // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
    const yMax = (histogramCache as any).bins.reduce((l: any, r: any) =>
      l > r ? l : r
    );

    // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
    (histogramCache as any).y = d3
      .scaleLinear()
      .domain([0, yMax])
      .range([topMargin + newHeight, topMargin]);

    return histogramCache;
  }

  // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types --- FIXME: disabled temporarily on migrate to TS.
  createQuery() {
    // @ts-expect-error ts-migrate(2339) FIXME: Property 'isObs' does not exist on type 'Readonly<... Remove this comment to see the full error message
    const { isObs, isGeneSetSummary, field, setGenes, annoMatrix } = this.props;
    const { schema } = annoMatrix;
    if (isObs) {
      return ["obs", field];
    }
    const varIndex = schema?.annotations?.var?.index;
    if (!varIndex) return null;

    if (isGeneSetSummary) {
      return [
        "X",
        {
          summarize: {
            method: "mean",
            field: "var",
            column: varIndex,
            values: [...setGenes.keys()],
          },
        },
      ];
    }

    // else, we assume it is a gene expression
    return [
      "X",
      {
        where: {
          field: "var",
          column: varIndex,
          value: field,
        },
      },
    ];
  }

  // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types --- FIXME: disabled temporarily on migrate to TS.
  render() {
    const {
      // @ts-expect-error ts-migrate(2339) FIXME: Property 'dispatch' does not exist on type 'Readon... Remove this comment to see the full error message
      dispatch,
      // @ts-expect-error ts-migrate(2339) FIXME: Property 'annoMatrix' does not exist on type 'Read... Remove this comment to see the full error message
      annoMatrix,
      // @ts-expect-error ts-migrate(2339) FIXME: Property 'field' does not exist on type 'Readonly<... Remove this comment to see the full error message
      field,
      // @ts-expect-error ts-migrate(2339) FIXME: Property 'isColorAccessor' does not exist on type ... Remove this comment to see the full error message
      isColorAccessor,
      // @ts-expect-error ts-migrate(2339) FIXME: Property 'isUserDefined' does not exist on type 'R... Remove this comment to see the full error message
      isUserDefined,
      // @ts-expect-error ts-migrate(2339) FIXME: Property 'isGeneSetSummary' does not exist on type... Remove this comment to see the full error message
      isGeneSetSummary,
      // @ts-expect-error ts-migrate(2339) FIXME: Property 'isScatterplotXXaccessor' does not exist ... Remove this comment to see the full error message
      isScatterplotXXaccessor,
      // @ts-expect-error ts-migrate(2339) FIXME: Property 'isScatterplotYYaccessor' does not exist ... Remove this comment to see the full error message
      isScatterplotYYaccessor,
      // @ts-expect-error ts-migrate(2339) FIXME: Property 'zebra' does not exist on type 'Readonly<... Remove this comment to see the full error message
      zebra,
      // @ts-expect-error ts-migrate(2339) FIXME: Property 'continuousSelectionRange' does not exist... Remove this comment to see the full error message
      continuousSelectionRange,
      // @ts-expect-error ts-migrate(2339) FIXME: Property 'isObs' does not exist on type 'Readonly<... Remove this comment to see the full error message
      isObs,
      // @ts-expect-error ts-migrate(2339) FIXME: Property 'mini' does not exist on type 'Readonly<{... Remove this comment to see the full error message
      mini,
      // @ts-expect-error ts-migrate(2339) FIXME: Property 'setGenes' does not exist on type 'Readon... Remove this comment to see the full error message
      setGenes,
    } = this.props;

    // @ts-expect-error ts-migrate(2339) FIXME: Property 'width' does not exist on type 'Readonly<... Remove this comment to see the full error message
    let { width } = this.props;
    if (!width) {
      width = mini ? WIDTH_MINI : WIDTH;
    }

    const fieldForId = field.replace(/\s/g, "_");
    const showScatterPlot = isUserDefined;

    let testClass = "histogram-continuous-metadata";
    if (isUserDefined) testClass = "histogram-user-gene";
    else if (isGeneSetSummary) testClass = "histogram-gene-set-summary";

    return (
      <Async
        watchFn={HistogramBrush.watchAsync}
        promiseFn={this.fetchAsyncProps}
        watchProps={{ annoMatrix, setGenes }}
      >
        <Async.Pending initial>
          <StillLoading displayName={field} zebra={zebra} />
        </Async.Pending>
        <Async.Rejected>
          {(error) => (
            <ErrorLoading zebra={zebra} error={error} displayName={field} />
          )}
        </Async.Rejected>
        <Async.Fulfilled>
          {(asyncProps) =>
            // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
            (asyncProps as any).OK2Render ? (
              <div
                id={`histogram_${fieldForId}`}
                data-testid={`histogram-${field}`}
                data-testclass={testClass}
                style={{
                  padding: mini ? 0 : globals.leftSidebarSectionPadding,
                  backgroundColor: zebra ? globals.lightestGrey : "white",
                }}
              >
                {!mini && isObs ? (
                  <HistogramHeader
                    // @ts-expect-error ts-migrate(2322) FIXME: Type '{ fieldId: any; isColorBy: any; isObs: any; ... Remove this comment to see the full error message
                    fieldId={field}
                    isColorBy={isColorAccessor}
                    isObs={isObs}
                    onColorByClick={this.handleColorAction(dispatch)}
                    onRemoveClick={isUserDefined ? this.removeHistogram : null}
                    isScatterPlotX={isScatterplotXXaccessor}
                    isScatterPlotY={isScatterplotYYaccessor}
                    onScatterPlotXClick={
                      showScatterPlot ? this.handleSetGeneAsScatterplotX : null
                    }
                    onScatterPlotYClick={
                      showScatterPlot ? this.handleSetGeneAsScatterplotY : null
                    }
                  />
                ) : null}
                <Histogram
                  field={field}
                  fieldForId={fieldForId}
                  // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
                  display={(asyncProps as any).isSingleValue ? "none" : "block"}
                  histogram={
                    mini
                      ? // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
                        (asyncProps as any).miniHistogram
                      : // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
                        (asyncProps as any).histogram
                  }
                  width={width}
                  height={mini ? HEIGHT_MINI : HEIGHT}
                  onBrush={this.onBrush}
                  onBrushEnd={this.onBrushEnd}
                  margin={mini ? MARGIN_MINI : MARGIN}
                  isColorBy={isColorAccessor}
                  selectionRange={continuousSelectionRange}
                  mini={mini}
                />
                {!mini && (
                  <HistogramFooter
                    // @ts-expect-error ts-migrate(2322) FIXME: Type '{ isGeneSetSummary: any; isObs: any; display... Remove this comment to see the full error message
                    isGeneSetSummary={isGeneSetSummary}
                    isObs={isObs}
                    displayName={field}
                    // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
                    hideRanges={(asyncProps as any).isSingleValue}
                    // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
                    rangeMin={(asyncProps as any).unclippedRange[0]}
                    // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
                    rangeMax={(asyncProps as any).unclippedRange[1]}
                    // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
                    rangeColorMin={(asyncProps as any).unclippedRangeColor[0]}
                    // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
                    rangeColorMax={(asyncProps as any).unclippedRangeColor[1]}
                  />
                )}
              </div>
            ) : null
          }
        </Async.Fulfilled>
      </Async>
    );
  }
}

export default HistogramBrush;
