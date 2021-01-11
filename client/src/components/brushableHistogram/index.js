/*
https://bl.ocks.org/mbostock/4341954
https://bl.ocks.org/mbostock/34f08d5e11952a80609169b7917d4172
https://bl.ocks.org/SpaceActuary/2f004899ea1b2bd78d6f1dbb2febf771
https://bl.ocks.org/mbostock/3019563
*/
import React, { useEffect, useRef, useState, useCallback } from "react";
import { Button, ButtonGroup, Icon, Tooltip } from "@blueprintjs/core";
import { connect } from "react-redux";
import * as d3 from "d3";
import { interpolateCool } from "d3-scale-chromatic";
import Async from "react-async";
import memoize from "memoize-one";
import { IconNames } from "@blueprintjs/icons";
import * as globals from "../../globals";
import actions from "../../actions";
import { histogramContinuous } from "../../util/dataframe/histogram";
import { makeContinuousDimensionName } from "../../util/nameCreators";
import significantDigits from "../../util/significantDigits";

function clamp(val, rng) {
  return Math.max(Math.min(val, rng[1]), rng[0]);
}

function maybeScientific(x) {
  let format = ",";
  const _ticks = x.ticks(4);

  if (x.domain().some((n) => Math.abs(n) >= 10000)) {
    /*
      heuristic: if the last tick d3 wants to render has one significant
      digit ie., 2000, render 2e+3, but if it's anything else ie., 42000000 render
      4.20e+n
    */
    format = significantDigits(_ticks[_ticks.length - 1]) === 1 ? ".0e" : ".2e";
  }

  return format;
}

const StillLoading = ({ zebra, displayName }) => {
  /*
  Render a loading indicator for the field.
  */
  return (
    <div
      style={{
        padding: globals.leftSidebarSectionPadding,
        backgroundColor: zebra ? globals.lightestGrey : "white",
      }}
    >
      <div
        style={{
          display: "flex",
          justifyContent: "space-between",
          justifyItems: "center",
          alignItems: "center",
        }}
      >
        <div style={{ minWidth: 30 }} />
        <div style={{ display: "flex", alignSelf: "center" }}>
          <span style={{ fontStyle: "italic" }}>{displayName}</span>
        </div>
        <div
          style={{
            display: "flex",
            justifyContent: "flex-end",
          }}
        >
          <Button minimal loading intent="primary" />
        </div>
      </div>
    </div>
  );
};

const ErrorLoading = ({ displayName, error, zebra }) => {
  console.log(error); // log to console as this is unexpected
  return (
    <div
      style={{
        padding: globals.leftSidebarSectionPadding,
        backgroundColor: zebra ? globals.lightestGrey : "white",
      }}
    >
      <span>{`Failure loading ${displayName}`}</span>
    </div>
  );
};

const HistogramFooter = React.memo(
  ({
    displayName,
    hideRanges,
    rangeMin,
    rangeMax,
    rangeColorMin,
    rangeColorMax,
    logFoldChange,
    pvalAdj,
  }) => {
    /*
  Footer of each histogram.  Will render range, title, and optionally
  differential expression info.

  Required props:
    * displayName - the displayName, aka "n_genes", "FOXP2", etc.
    * hideRanges - true/false, enables/disable rendering of ranges
    * range - length two array, [min, max], containing the range values to display
    * rangeColor - length two array, [mincolor, maxcolor], each a CSS color
    * logFoldChange - lfc to display, optional.
    * pValue - pValue to display, optional.
  */
    return (
      <div>
        <div
          style={{
            display: "flex",
            justifyContent: hideRanges ? "center" : "space-between",
          }}
        >
          <span
            style={{
              color: rangeColorMin,
              display: hideRanges ? "none" : "block",
            }}
          >
            min {rangeMin.toPrecision(4)}
          </span>
          <span
            data-testclass="brushable-histogram-field-name"
            style={{ fontStyle: "italic" }}
          >
            {displayName}
          </span>
          <div style={{ display: hideRanges ? "block" : "none" }}>
            : {rangeMin}
          </div>
          <span
            style={{
              color: rangeColorMax,
              display: hideRanges ? "none" : "block",
            }}
          >
            max {rangeMax.toPrecision(4)}
          </span>
        </div>

        {logFoldChange !== undefined && pvalAdj !== undefined ? (
          <div
            style={{
              display: "flex",
              justifyContent: "center",
              alignItems: "baseline",
            }}
          >
            <span>
              <strong>log fold change:</strong>
              {` ${logFoldChange.toPrecision(4)}`}
            </span>
            <span
              style={{
                marginLeft: 7,
                padding: 2,
              }}
            >
              <strong>p-value (adj):</strong>
              {pvalAdj < 0.0001 ? " < 0.0001" : ` ${pvalAdj.toFixed(4)}`}
            </span>
          </div>
        ) : null}
      </div>
    );
  }
);

const HistogramHeader = React.memo(
  ({
    fieldId,
    isColorBy,
    onColorByClick,
    onRemoveClick,
    isScatterPlotX,
    isScatterPlotY,
    onScatterPlotXClick,
    onScatterPlotYClick,
    isObs,
  }) => {
    /*
      Render the toolbar for the histogram.  Props:
        * fieldId - field identifier, used for various IDs
        * isColorBy - true/false, is this the current color-by
        * onColorByClick - color-by click handler
        * onRemoveClick - optional handler for remove.  Button will not render if not defined.
        * isScatterPlotX - optional, true/false if currently the X scatterplot field
        * isScatterPlotY - optional, true/false if currently the Y scatterplot field
        * onScatterPlotXClick - optional, handler for scatterPlot X button.
        * onScatterPlotYClick - optional, handler for scatterPlot X button.

      Scatterplot controls will not render if either handler unspecified.
    */

    const memoizedColorByCallback = useCallback(
      () => onColorByClick(fieldId, isObs),
      [fieldId, isObs]
    );

    return (
      <div
        style={{
          display: "flex",
          justifyContent: "flex-end",
          paddingBottom: "8px",
        }}
      >
        {onScatterPlotXClick && onScatterPlotYClick ? (
          <span>
            <Icon icon={IconNames.SCATTER_PLOT} style={{ marginRight: 7 }} />
            <ButtonGroup style={{ marginRight: 7 }}>
              <Button
                data-testid={`plot-x-${fieldId}`}
                onClick={onScatterPlotXClick}
                active={isScatterPlotX}
                intent={isScatterPlotX ? "primary" : "none"}
              >
                plot x
              </Button>
              <Button
                data-testid={`plot-y-${fieldId}`}
                onClick={onScatterPlotYClick}
                active={isScatterPlotY}
                intent={isScatterPlotY ? "primary" : "none"}
              >
                plot y
              </Button>
            </ButtonGroup>
          </span>
        ) : null}
        {onRemoveClick ? (
          <Button
            minimal
            onClick={onRemoveClick}
            style={{
              color: globals.blue,
              cursor: "pointer",
              marginLeft: 7,
            }}
          >
            remove
          </Button>
        ) : null}
        <Tooltip
          content="Use as color scale"
          position="bottom"
          hoverOpenDelay={globals.tooltipHoverOpenDelay}
        >
          <Button
            onClick={memoizedColorByCallback}
            active={isColorBy}
            intent={isColorBy ? "primary" : "none"}
            data-testclass="colorby"
            data-testid={`colorby-${fieldId}`}
            icon="tint"
          />
        </Tooltip>
      </div>
    );
  }
);

const Histogram = ({
  field,
  fieldForId,
  display,
  histogram,
  width,
  height,
  onBrush,
  onBrushEnd,
  margin,
  isColorBy,
  selectionRange,
}) => {
  const svgRef = useRef(null);
  const [brush, setBrush] = useState(null);

  useEffect(() => {
    /*
    Create the d3 histogram
    */
    const { marginLeft, marginRight, marginBottom, marginTop } = margin;
    const { x, y, bins, binStart, binEnd, binWidth } = histogram;
    const svg = d3.select(svgRef.current);

    /* Remove everything */
    svg.selectAll("*").remove();

    /* Set margins within the SVG */
    const container = svg
      .attr("width", width + marginLeft + marginRight)
      .attr("height", height + marginTop + marginBottom)
      .append("g")
      .attr("class", "histogram-container")
      .attr("transform", `translate(${marginLeft},${marginTop})`);

    const colorScale = d3
      .scaleSequential(interpolateCool)
      .domain([0, bins.length]);

    const histogramScale = d3
      .scaleLinear()
      .domain(x.domain())
      .range([
        colorScale.domain()[1],
        colorScale.domain()[0],
      ]); /* we flip this to make colors dark if high in the color scale */

    if (binWidth > 0) {
      /* BINS */
      container
        .insert("g", "*")
        .selectAll("rect")
        .data(bins)
        .enter()
        .append("rect")
        .attr("x", (d, i) => x(binStart(i)) + 1)
        .attr("y", (d) => y(d))
        .attr("width", (d, i) => x(binEnd(i)) - x(binStart(i)) - 1)
        .attr("height", (d) => y(0) - y(d))
        .style(
          "fill",
          isColorBy ? (d, i) => colorScale(histogramScale(binStart(i))) : "#bbb"
        );
    }

    // BRUSH
    // Note the brushable area is bounded by the data on three sides, but goes down to cover the x-axis
    const brushX = d3
      .brushX()
      .extent([
        [x.range()[0], y.range()[1]],
        [x.range()[1], marginTop + height + marginBottom],
      ])
      /*
      emit start so that the Undoable history can save an undo point
      upon drag start, and ignore the subsequent intermediate drag events.
      */
      .on("start", onBrush(field, x.invert, "start"))
      .on("brush", onBrush(field, x.invert, "brush"))
      .on("end", onBrushEnd(field, x.invert));

    const brushXselection = container
      .insert("g")
      .attr("class", "brush")
      .attr("data-testid", `${svgRef.current.dataset.testid}-brushable-area`)
      .call(brushX);

    /* X AXIS */
    container
      .insert("g")
      .attr("class", "axis axis--x")
      .attr("transform", `translate(0,${marginTop + height})`)
      .call(
        d3
          .axisBottom(x)
          .ticks(4)
          .tickFormat(d3.format(maybeScientific(x)))
      );

    /* Y AXIS */
    container
      .insert("g")
      .attr("class", "axis axis--y")
      .attr("transform", `translate(${marginLeft + width},0)`)
      .call(
        d3
          .axisRight(y)
          .ticks(3)
          .tickFormat(
            d3.format(
              y.domain().some((n) => Math.abs(n) >= 10000) ? ".0e" : ","
            )
          )
      );

    /* axis style */
    svg.selectAll(".axis text").style("fill", "rgb(80,80,80)");
    svg.selectAll(".axis path").style("stroke", "rgb(230,230,230)");
    svg.selectAll(".axis line").style("stroke", "rgb(230,230,230)");

    setBrush({ brushX, brushXselection });
  }, [histogram, isColorBy]);

  useEffect(() => {
    /*
    paint/update selection brush
    */
    if (!brush) return;
    const { brushX, brushXselection } = brush;
    const selection = d3.brushSelection(brushXselection.node());
    if (!selectionRange && selection) {
      /* no active selection - clear brush */
      brushXselection.call(brushX.move, null);
    } else if (selectionRange) {
      const { x, domain } = histogram;
      const [min, max] = domain;
      const x0 = x(clamp(selectionRange[0], [min, max]));
      const x1 = x(clamp(selectionRange[1], [min, max]));
      if (!selection) {
        /* there is an active selection, but no brush - set the brush */
        brushXselection.call(brushX.move, [x0, x1]);
      } else {
        /* there is an active selection and a brush - make sure they match */
        const moveDeltaThreshold = 1;
        const dX0 = Math.abs(x0 - selection[0]);
        const dX1 = Math.abs(x1 - selection[1]);
        /*
        only update the brush if it is grossly incorrect,
        as defined by the moveDeltaThreshold
        */
        if (dX0 > moveDeltaThreshold || dX1 > moveDeltaThreshold) {
          brushXselection.call(brushX.move, [x0, x1]);
        }
      }
    }
  }, [brush, selectionRange]);

  return (
    <svg
      style={{ display }}
      width={width}
      height={height}
      id={`histogram_${fieldForId}_svg`}
      data-testclass="histogram-plot"
      data-testid={`histogram-${field}-plot`}
      ref={svgRef}
    />
  );
};

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
    this.margin = {
      marginLeft,
      marginRight,
      marginBottom,
      marginTop,
    };
    this.width = 340 - marginLeft - marginRight;
    this.height = 135 - marginTop - marginBottom;
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

    const histogram = this.calcHistogramCache(
      column,
      this.margin,
      this.width,
      this.height
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
      range,
      unclippedRange,
      unclippedRangeColor,
      isSingleValue,
      OK2Render,
    };
  };

  // eslint-disable-next-line class-methods-use-this -- instance method allows for memoization per annotation
  calcHistogramCache(col, margin, width, height) {
    /*
     recalculate expensive stuff, notably bins, summaries, etc.
    */
    const histogramCache = {};
    const summary = col.summarize();
    const { min: domainMin, max: domainMax } = summary;
    const numBins = 40;
    const { marginTop, marginLeft } = margin;

    histogramCache.domain = [domainMin, domainMax];

    histogramCache.x = d3
      .scaleLinear()
      .domain([domainMin, domainMax])
      .range([marginLeft, marginLeft + width]);

    histogramCache.bins = histogramContinuous(col, numBins, [
      domainMin,
      domainMax,
    ]);
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
    } = this.props;
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
                  padding: globals.leftSidebarSectionPadding,
                  backgroundColor: zebra ? globals.lightestGrey : "white",
                }}
              >
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
                <Histogram
                  field={field}
                  fieldForId={fieldForId}
                  display={asyncProps.isSingleValue ? "none" : "block"}
                  histogram={asyncProps.histogram}
                  width={this.width}
                  height={this.height}
                  onBrush={this.onBrush}
                  onBrushEnd={this.onBrushEnd}
                  margin={this.margin}
                  isColorBy={isColorAccessor}
                  selectionRange={continuousSelectionRange}
                />
                <HistogramFooter
                  displayName={field}
                  hideRanges={asyncProps.isSingleValue}
                  rangeMin={asyncProps.unclippedRange[0]}
                  rangeMax={asyncProps.unclippedRange[1]}
                  rangeColorMin={asyncProps.unclippedRangeColor[0]}
                  rangeColorMax={asyncProps.unclippedRangeColor[1]}
                  logFoldChange={logFoldChange}
                  pvalAdj={pvalAdj}
                />
              </div>
            ) : null
          }
        </Async.Fulfilled>
      </Async>
    );
  }
}

export default HistogramBrush;
