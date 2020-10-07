/*
https://bl.ocks.org/mbostock/4341954
https://bl.ocks.org/mbostock/34f08d5e11952a80609169b7917d4172
https://bl.ocks.org/SpaceActuary/2f004899ea1b2bd78d6f1dbb2febf771
https://bl.ocks.org/mbostock/3019563
*/
import React from "react";
import { connect } from "react-redux";
import * as d3 from "d3";
import Async from "react-async";
import memoize from "memoize-one";
import * as globals from "../../globals";
import actions from "../../actions";
import { histogramContinuous } from "../../util/dataframe/histogram";
import { makeContinuousDimensionName } from "../../util/nameCreators";
import HistogramHeader from "./header";
import Histogram from "./histogram";
import HistogramFooter from "./footer";
import StillLoading from "./loading";
import ErrorLoading from "./error";

@connect((state, ownProps) => {
  const { isObs, isUserDefined, isDiffExp, field } = ownProps;
  const myName = makeContinuousDimensionName(
    { isObs, isUserDefined, isDiffExp },
    field
  );
  return {
    annoMatrix: state.annoMatrix,
    isScatterplotXXaccessor: state.controls.scatterplotXXaccessor === field,
    isScatterplotYYaccessor: state.controls.scatterplotYYaccessor === field,
    continuousSelectionRange: state.continuousSelection[myName],
    isColorAccessor: state.colors.colorAccessor === field,
  };
})
class HistogramBrush extends React.PureComponent {
  /* memoized closure to prevent HistogramHeader unecessary repaint */
  handleColorAction = memoize((dispatch) => (field, isObs) => {
    if (isObs) {
      dispatch({
        type: "color by continuous metadata",
        colorAccessor: field,
      });
    } else {
      dispatch(actions.requestSingleGeneExpressionCountsForColoringPOST(field));
    }
  });

  constructor(props) {
    super(props);

    const marginLeft = 10; // Space for 0 tick label on X axis
    const marginRight = 54; // space for Y axis & labels
    const marginBottom = 25; // space for X axis & labels
    const marginTop = 3;

    this.state = {
      margin: {
        marginLeft,
        marginRight,
        marginBottom,
        marginTop,
      },
      width: 340 - marginLeft - marginRight,
      height: 135 - marginTop - marginBottom,
      marginMini: {
        marginLeft: 0, // Space for 0 tick label on X axis
        marginRight: 0, // space for Y axis & labels
        marginBottom: 0, // space for X axis & labels
        marginTop: 0,
      },
      widthMini: 120,
      heightMini: 15,
    };
  }

  onBrush = (selection, x, eventType) => {
    const type = `continuous metadata histogram ${eventType}`;
    return () => {
      const { dispatch, field, isObs, isUserDefined, isDiffExp } = this.props;

      // ignore programmatically generated events
      if (!d3.event.sourceEvent) return;
      // ignore cascading events, which are programmatically generated
      if (d3.event.sourceEvent.sourceEvent) return;

      const query = this.createQuery();
      const range = d3.event.selection
        ? [x(d3.event.selection[0]), x(d3.event.selection[1])]
        : null;
      const otherProps = {
        selection: field,
        continuousNamespace: {
          isObs,
          isUserDefined,
          isDiffExp,
        },
      };
      dispatch(
        actions.selectContinuousMetadataAction(type, query, range, otherProps)
      );
    };
  };

  onBrushEnd = (selection, x) => {
    return () => {
      const { dispatch, field, isObs, isUserDefined, isDiffExp } = this.props;
      const minAllowedBrushSize = 10;
      const smallAmountToAvoidInfiniteLoop = 0.1;

      // ignore programmatically generated events
      if (!d3.event.sourceEvent) return;
      // ignore cascading events, which are programmatically generated
      if (d3.event.sourceEvent.sourceEvent) return;

      let type;
      let range = null;
      if (d3.event.selection) {
        type = "continuous metadata histogram end";
        if (
          d3.event.selection[1] - d3.event.selection[0] >
          minAllowedBrushSize
        ) {
          range = [x(d3.event.selection[0]), x(d3.event.selection[1])];
        } else {
          /* the user selected range is too small and will be hidden #587, so take control of it procedurally */
          /* https://stackoverflow.com/questions/12354729/d3-js-limit-size-of-brush */

          const procedurallyResizedBrushWidth =
            d3.event.selection[0] +
            minAllowedBrushSize +
            smallAmountToAvoidInfiniteLoop; //

          range = [x(d3.event.selection[0]), x(procedurallyResizedBrushWidth)];
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
          isDiffExp,
        },
      };
      dispatch(
        actions.selectContinuousMetadataAction(type, query, range, otherProps)
      );
    };
  };

  handleSetGeneAsScatterplotX = () => {
    const { dispatch, field } = this.props;
    dispatch({
      type: "set scatterplot x",
      data: field,
    });
  };

  handleSetGeneAsScatterplotY = () => {
    const { dispatch, field } = this.props;
    dispatch({
      type: "set scatterplot y",
      data: field,
    });
  };

  removeHistogram = () => {
    const {
      dispatch,
      field,
      isColorAccessor,
      isScatterplotXXaccessor,
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

  fetchAsyncProps = async () => {
    const { annoMatrix } = this.props;
    const {
      margin,
      width,
      height,
      marginMini,
      widthMini,
      heightMini,
    } = this.state;
    const { isClipped } = annoMatrix;

    const query = this.createQuery();
    const df = await annoMatrix.fetch(...query);
    const column = df.icol(0);

    // if we are clipped, fetch both our value and our unclipped value,
    // as we need the absolute min/max range, not just the clipped min/max.
    const summary = column.summarize();
    const range = [summary.min, summary.max];

    let unclippedRange = [...range];
    if (isClipped) {
      const parent = await annoMatrix.viewOf.fetch(...query);
      const { min, max } = parent.icol(0).summarize();
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

    const histogram = this.calcHistogramCache(column, margin, width, height);
    const miniHistogram = this.calcHistogramCache(
      column,
      marginMini,
      widthMini,
      heightMini
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

  // eslint-disable-next-line class-methods-use-this -- instance method allows for memoization per annotation
  calcHistogramCache(col, margin, width, height) {
    /* make this more structured and doing a forEach */
    /*
     recalculate expensive stuff, notably bins, summaries, etc.
    */
    const histogramCache = {}; /* maybe change this so that it computes ... */
    const summary = col.summarize(); /* this is memoized, so it's free the second time you call it */
    const { min: domainMin, max: domainMax } = summary;
    const numBins = 40;
    const { marginTop, marginLeft } = margin; /* changes with mini */

    histogramCache.domain = [
      domainMin,
      domainMax,
    ]; /* doesn't change with mini */

    histogramCache.x = d3
      .scaleLinear()
      .domain([domainMin, domainMax])
      .range([marginLeft, marginLeft + width]);

    histogramCache.bins = histogramContinuous(col, numBins, [
      domainMin,
      domainMax,
    ]); /* memoized */
    histogramCache.binWidth = (domainMax - domainMin) / numBins;

    histogramCache.binStart = (i) => domainMin + i * histogramCache.binWidth;
    histogramCache.binEnd = (i) =>
      domainMin + (i + 1) * histogramCache.binWidth;

    const yMax = histogramCache.bins.reduce((l, r) => (l > r ? l : r));

    histogramCache.y = d3
      .scaleLinear()
      .domain([0, yMax])
      .range([marginTop + height, marginTop]);

    return histogramCache;
  }

  createQuery() {
    const { isObs, field, annoMatrix } = this.props;
    const { schema } = annoMatrix;
    if (isObs) {
      return ["obs", field];
    }
    const varIndex = schema?.annotations?.var?.index;
    if (!varIndex) return null;
    return [
      "X",
      {
        field: "var",
        column: varIndex,
        value: field,
      },
    ];
  }

  render() {
    const {
      dispatch,
      annoMatrix,
      field,
      isColorAccessor,
      isUserDefined,
      isDiffExp,
      logFoldChange,
      pvalAdj,
      isScatterplotXXaccessor,
      isScatterplotYYaccessor,
      zebra,
      continuousSelectionRange,
      isObs,
      mini,
    } = this.props;
    const {
      margin,
      width,
      height,
      marginMini,
      widthMini,
      heightMini,
    } = this.state;
    const fieldForId = field.replace(/\s/g, "_");
    const showScatterPlot = isDiffExp || isUserDefined;

    return (
      <Async watch={annoMatrix} promiseFn={this.fetchAsyncProps}>
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
            asyncProps.OK2Render ? (
              <div
                id={`histogram_${fieldForId}`}
                data-testid={`histogram-${field}`}
                data-testclass={
                  isDiffExp
                    ? "histogram-diffexp"
                    : isUserDefined
                    ? "histogram-user-gene"
                    : "histogram-continuous-metadata"
                }
                style={{
                  padding: mini ? 0 : globals.leftSidebarSectionPadding,
                  backgroundColor: zebra ? globals.lightestGrey : "white",
                }}
              >
                {!mini && isObs ? (
                  <HistogramHeader
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
                  display={asyncProps.isSingleValue ? "none" : "block"}
                  histogram={
                    mini ? asyncProps.miniHistogram : asyncProps.histogram
                  }
                  width={mini ? widthMini : width}
                  height={mini ? heightMini : height}
                  onBrush={this.onBrush}
                  onBrushEnd={this.onBrushEnd}
                  margin={mini ? marginMini : margin}
                  isColorBy={isColorAccessor}
                  selectionRange={continuousSelectionRange}
                  mini={mini}
                />
                {!mini ? (
                  <HistogramFooter
                    isObs={isObs}
                    displayName={field}
                    hideRanges={asyncProps.isSingleValue}
                    rangeMin={asyncProps.unclippedRange[0]}
                    rangeMax={asyncProps.unclippedRange[1]}
                    rangeColorMin={asyncProps.unclippedRangeColor[0]}
                    rangeColorMax={asyncProps.unclippedRangeColor[1]}
                    logFoldChange={logFoldChange}
                    pvalAdj={pvalAdj}
                  />
                ) : null}
              </div>
            ) : null
          }
        </Async.Fulfilled>
      </Async>
    );
  }
}

export default HistogramBrush;
