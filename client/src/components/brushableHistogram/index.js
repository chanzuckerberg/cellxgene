/*
https://bl.ocks.org/mbostock/4341954
https://bl.ocks.org/mbostock/34f08d5e11952a80609169b7917d4172
https://bl.ocks.org/SpaceActuary/2f004899ea1b2bd78d6f1dbb2febf771
https://bl.ocks.org/mbostock/3019563
*/
// jshint esversion: 6
import React from "react";
import { Button, ButtonGroup, Tooltip } from "@blueprintjs/core";
import { connect } from "react-redux";
import * as d3 from "d3";
import memoize from "memoize-one";
import { interpolateCool } from "d3-scale-chromatic";
import * as globals from "../../globals";
import actions from "../../actions";
import { histogramContinuous } from "../../util/dataframe/histogram";
import { makeContinuousDimensionName } from "../../util/nameCreators";

function clamp(val, rng) {
  return Math.max(Math.min(val, rng[1]), rng[0]);
}

@connect((state, ownProps) => {
  const { isObs, isUserDefined, isDiffExp, field } = ownProps;
  const myName = makeContinuousDimensionName(
    { isObs, isUserDefined, isDiffExp },
    field
  );
  return {
    world: state.world,
    isScatterplotXXaccessor: state.controls.scatterplotXXaccessor === field,
    isScatterplotYYaccessor: state.controls.scatterplotYYaccessor === field,
    continuousSelectionRange: state.continuousSelection[myName],
    isColorAccessor: state.colors.colorAccessor === field,
  };
})
class HistogramBrush extends React.PureComponent {
  static getColumn(world, field, clipped = true) {
    /*
    Return the underlying Dataframe column for our field.   By default,
    returns the clipped column.   If clipped===false, will return the
    unclipped column.
    */
    const obsAnnotations = clipped
      ? world.obsAnnotations
      : world.unclipped.obsAnnotations;
    const varData = clipped ? world.varData : world.unclipped.varData;
    if (obsAnnotations.hasCol(field)) {
      return obsAnnotations.col(field);
    }
    return varData.col(field);
  }

  calcHistogramCache = memoize((col) => {
    /*
     recalculate expensive stuff, notably bins, summaries, etc.
    */
    const histogramCache = {};
    const summary = col.summarize();
    const { min: domainMin, max: domainMax } = summary;
    const numBins = 40;

    histogramCache.x = d3
      .scaleLinear()
      .domain([domainMin, domainMax])
      .range([this.marginLeft, this.marginLeft + this.width]);

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
      .range([this.marginTop + this.height, this.marginTop]);

    return histogramCache;
  });

  constructor(props) {
    super(props);

    this.marginLeft = 10; // Space for 0 tick label on X axis
    this.marginRight = 54; // space for Y axis & labels
    this.marginBottom = 25; // space for X axis & labels
    this.marginTop = 3;

    this.width = 340 - this.marginLeft - this.marginRight;
    this.height = 135 - this.marginTop - this.marginBottom;
  }

  componentDidMount() {
    const { field, isColorAccessor } = this.props;

    this.renderHistogram(this._histogram, field, isColorAccessor);
  }

  componentDidUpdate(prevProps) {
    const {
      field,
      world,
      continuousSelectionRange: range,
      isColorAccessor,
    } = this.props;
    const { x } = this._histogram;
    let { brushX, brushXselection } = this.state;

    const dfColumn = HistogramBrush.getColumn(world, field);
    const oldDfColumn = HistogramBrush.getColumn(
      prevProps.world,
      prevProps.field
    );

    const rangeChanged = range !== prevProps.continuousSelectionRange;
    const dfChanged = dfColumn !== oldDfColumn;
    const colorSelectionChanged = prevProps.isColorAccessor !== isColorAccessor;

    if (dfChanged || colorSelectionChanged) {
      ({ brushX, brushXselection } = this.renderHistogram(
        this._histogram,
        field,
        isColorAccessor
      ));
    }

    /*
    if the selection has changed, ensure that the brush correctly reflects
    the underlying selection.
    */
    if (
      (dfChanged || rangeChanged || colorSelectionChanged) &&
      brushXselection
    ) {
      const selection = d3.brushSelection(brushXselection.node());
      if (!range && selection) {
        /* no active selection - clear brush */
        brushXselection.call(brushX.move, null);
      } else if (range) {
        const { min, max } = dfColumn.summarize();
        const x0 = x(clamp(range[0], [min, max]));
        const x1 = x(clamp(range[1], [min, max]));
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
    }
  }

  onBrush(selection, x, eventType) {
    const type = `continuous metadata histogram ${eventType}`;
    return () => {
      const { dispatch, field, isObs, isUserDefined, isDiffExp } = this.props;

      // ignore programmatically generated events
      if (!d3.event.sourceEvent) return;
      // ignore cascading events, which are programmatically generated
      if (d3.event.sourceEvent.sourceEvent) return;

      if (d3.event.selection) {
        dispatch({
          type,
          selection: field,
          continuousNamespace: {
            isObs,
            isUserDefined,
            isDiffExp,
          },
          range: [x(d3.event.selection[0]), x(d3.event.selection[1])],
        });
      } else {
        dispatch({
          type,
          selection: field,
          continuousNamespace: {
            isObs,
            isUserDefined,
            isDiffExp,
          },
          range: null,
        });
      }
    };
  }

  onBrushEnd(selection, x) {
    return () => {
      const { dispatch, field, isObs, isUserDefined, isDiffExp } = this.props;
      const minAllowedBrushSize = 10;
      const smallAmountToAvoidInfiniteLoop = 0.1;

      // ignore programmatically generated events
      if (!d3.event.sourceEvent) return;
      // ignore cascading events, which are programmatically generated
      if (d3.event.sourceEvent.sourceEvent) return;

      if (d3.event.selection) {
        let _range;

        if (
          d3.event.selection[1] - d3.event.selection[0] >
          minAllowedBrushSize
        ) {
          _range = [x(d3.event.selection[0]), x(d3.event.selection[1])];
        } else {
          /* the user selected range is too small and will be hidden #587, so take control of it procedurally */
          /* https://stackoverflow.com/questions/12354729/d3-js-limit-size-of-brush */

          const procedurallyResizedBrushWidth =
            d3.event.selection[0] +
            minAllowedBrushSize +
            smallAmountToAvoidInfiniteLoop; //

          _range = [x(d3.event.selection[0]), x(procedurallyResizedBrushWidth)];
        }

        dispatch({
          type: "continuous metadata histogram end",
          selection: field,
          continuousNamespace: {
            isObs,
            isUserDefined,
            isDiffExp,
          },
          range: _range,
        });
      } else {
        dispatch({
          type: "continuous metadata histogram cancel",
          selection: field,
          continuousNamespace: {
            isObs,
            isUserDefined,
            isDiffExp,
          },
        });
      }
    };
  }

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

  handleColorAction = () => {
    const { dispatch, field, world, ranges } = this.props;

    if (world.obsAnnotations.hasCol(field)) {
      dispatch({
        type: "color by continuous metadata",
        colorAccessor: field,
        rangeForColorAccessor: ranges,
      });
    } else if (world.varData.hasCol(field)) {
      dispatch(actions.requestSingleGeneExpressionCountsForColoringPOST(field));
    }
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

  drawHistogram(svgRef) {
    const { field, world } = this.props;
    const col = HistogramBrush.getColumn(world, field);
    this._histogram = {
      ...this.calcHistogramCache(col),
      svgRef,
    };
  }

  renderHistogram(histogram, field, isColorAccessor) {
    const { x, y, bins, svgRef, binStart, binEnd, binWidth } = histogram;
    const svg = d3.select(svgRef);

    /* Remove everything */
    svg.selectAll("*").remove();

    /* Set margins within the SVG */
    const container = svg
      .attr("width", this.width + this.marginLeft + this.marginRight)
      .attr("height", this.height + this.marginTop + this.marginBottom)
      .append("g")
      .attr("class", "histogram-container")
      .attr("transform", `translate(${this.marginLeft},${this.marginTop})`);

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
          isColorAccessor
            ? (d, i) => colorScale(histogramScale(binStart(i)))
            : "#bbb"
        );
    }

    // BRUSH
    // Note the brushable area is bounded by the data on three sides, but goes down to cover the x-axis
    const brushX = d3
      .brushX()
      .extent([
        [x.range()[0], y.range()[1]],
        [x.range()[1], this.marginTop + this.height + this.marginBottom],
      ])
      /*
      emit start so that the Undoable history can save an undo point
      upon drag start, and ignore the subsequent intermediate drag events.
      */
      .on("start", this.onBrush(field, x.invert, "start").bind(this))
      .on("brush", this.onBrush(field, x.invert, "brush").bind(this))
      .on("end", this.onBrushEnd(field, x.invert).bind(this));

    const brushXselection = container
      .insert("g")
      .attr("class", "brush")
      .attr("data-testid", `${svgRef.dataset.testid}-brushable-area`)
      .call(brushX);

    /* X AXIS */
    container
      .insert("g")
      .attr("class", "axis axis--x")
      .attr("transform", `translate(0,${this.marginTop + this.height})`)
      .call(
        d3
          .axisBottom(x)
          .ticks(4)
          .tickFormat(
            d3.format(
              x.domain().some((n) => Math.abs(n) >= 10000) ? ".2e" : ","
            )
          )
      );

    /* Y AXIS */
    container
      .insert("g")
      .attr("class", "axis axis--y")
      .attr("transform", `translate(${this.marginLeft + this.width},0)`)
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

    const newState = { brushX, brushXselection };
    this.setState(newState);
    return newState;
  }

  render() {
    const {
      field,
      world,
      isColorAccessor,
      isUserDefined,
      isDiffExp,
      logFoldChange,
      pvalAdj,
      isScatterplotXXaccessor,
      isScatterplotYYaccessor,
      zebra,
    } = this.props;
    const fieldForId = field.replace(/\s/g, "_");
    const {
      min: unclippedRangeMin,
      max: unclippedRangeMax,
    } = HistogramBrush.getColumn(world, field, false).summarize();
    const unclippedRangeMinColor =
      world.clipQuantiles.min === 0 ? "#bbb" : globals.blue;
    const unclippedRangeMaxColor =
      world.clipQuantiles.max === 1 ? "#bbb" : globals.blue;

    return (
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
        <div
          style={{
            display: "flex",
            justifyContent: "flex-end",
            paddingBottom: "8px",
          }}
        >
          {isDiffExp || isUserDefined ? (
            <span>
              <span
                style={{ marginRight: 7 }}
                className="bp3-icon-standard bp3-icon-scatter-plot"
              />
              <ButtonGroup style={{ marginRight: 7 }}>
                <Button
                  data-testid={`plot-x-${field}`}
                  onClick={this.handleSetGeneAsScatterplotX}
                  active={isScatterplotXXaccessor}
                  intent={isScatterplotXXaccessor ? "primary" : "none"}
                >
                  plot x
                </Button>
                <Button
                  data-testid={`plot-y-${field}`}
                  onClick={this.handleSetGeneAsScatterplotY}
                  active={isScatterplotYYaccessor}
                  intent={isScatterplotYYaccessor ? "primary" : "none"}
                >
                  plot y
                </Button>
              </ButtonGroup>
            </span>
          ) : null}
          {isUserDefined ? (
            <Button
              minimal
              onClick={this.removeHistogram}
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
              onClick={this.handleColorAction}
              active={isColorAccessor}
              intent={isColorAccessor ? "primary" : "none"}
              data-testclass="colorby"
              data-testid={`colorby-${field}`}
              icon="tint"
            />
          </Tooltip>
        </div>
        <svg
          width={this.width}
          height={this.height}
          id={`histogram_${fieldForId}_svg`}
          data-testclass="histogram-plot"
          data-testid={`histogram-${field}-plot`}
          ref={(svgRef) => this.drawHistogram(svgRef)}
        />
        <div
          style={{
            display: "flex",
            justifyContent: "space-between",
          }}
        >
          <span style={{ color: unclippedRangeMinColor }}>
            min {unclippedRangeMin.toPrecision(4)}
          </span>
          <span
            data-testclass="brushable-histogram-field-name"
            style={{ fontStyle: "italic" }}
          >
            {field}
          </span>
          <span style={{ color: unclippedRangeMaxColor }}>
            max {unclippedRangeMax.toPrecision(4)}
          </span>
        </div>

        {isDiffExp ? (
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
}

export default HistogramBrush;
