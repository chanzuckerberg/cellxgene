/*
https://bl.ocks.org/mbostock/4341954
https://bl.ocks.org/mbostock/34f08d5e11952a80609169b7917d4172
https://bl.ocks.org/SpaceActuary/2f004899ea1b2bd78d6f1dbb2febf771
*/
// jshint esversion: 6
import React from "react";
import _ from "lodash";
import { Button, ButtonGroup, Tooltip } from "@blueprintjs/core";
import { connect } from "react-redux";
import * as d3 from "d3";
import memoize from "memoize-one";
import * as globals from "../../globals";
import actions from "../../actions";
import finiteExtent from "../../util/finiteExtent";
import { makeContinuousDimensionName } from "../../util/nameCreators";

@connect(state => ({
  world: state.world,
  scatterplotXXaccessor: state.controls.scatterplotXXaccessor,
  scatterplotYYaccessor: state.controls.scatterplotYYaccessor,
  continuousSelection: state.continuousSelection,
  differential: state.differential,
  colorAccessor: state.colors.colorAccessor,
  obsAnnotations: _.get(state.world, "obsAnnotations", null)
}))
class HistogramBrush extends React.Component {
  calcHistogramCache = memoize((obsAnnotations, field, rangeMin, rangeMax) => {
    const { world } = this.props;
    const histogramCache = {};

    histogramCache.y = d3
      .scaleLinear()
      .range([this.height - this.marginBottom, 0]);

    if (obsAnnotations.hasCol(field)) {
      // recalculate expensive stuff
      const allValuesForContinuousFieldAsArray = obsAnnotations
        .col(field)
        .asArray();

      histogramCache.x = d3
        .scaleLinear()
        .domain([rangeMin, rangeMax])
        .range([0, this.width]);

      histogramCache.bins = d3
        .histogram()
        .domain(histogramCache.x.domain())
        .thresholds(40)(allValuesForContinuousFieldAsArray);

      histogramCache.numValues = allValuesForContinuousFieldAsArray.length;
    } else if (world.varData.hasCol(field)) {
      const varValues = world.varData.col(field).asArray();

      histogramCache.x = d3
        .scaleLinear()
        .domain(
          finiteExtent(varValues)
        ) /* replace this if we have ranges for genes back from server like we do for annotations on cells */
        .range([0, this.width]);

      histogramCache.bins = d3
        .histogram()
        .domain(histogramCache.x.domain())
        .thresholds(40)(varValues);

      histogramCache.numValues = varValues.length;
    }

    return histogramCache;
  });

  constructor(props) {
    super(props);

    this.width = 340;
    this.height = 100;
    this.marginBottom = 20;
  }

  componentDidMount() {
    const { field } = this.props;
    const { x, y, bins, numValues, svgRef } = this._histogram;

    this.renderAxesBrushBins(x, y, bins, numValues, svgRef, field);
  }

  componentDidUpdate(prevProps) {
    const { field, obsAnnotations, continuousSelection } = this.props;
    const { x, y, bins, numValues, svgRef } = this._histogram;

    if (obsAnnotations !== prevProps.obsAnnotations) {
      this.renderAxesBrushBins(x, y, bins, numValues, svgRef, field);
    }

    /*
    if the selection has changed, ensure that the brush correctly reflects
    the underlying selection.
    */
    if (continuousSelection !== prevProps.continuousSelection) {
      const { isObs, isUserDefined, isDiffExp } = this.props;
      const myName = makeContinuousDimensionName(
        { isObs, isUserDefined, isDiffExp },
        field
      );
      const range = continuousSelection[myName];
      const { brushXselection, brushX } = this.state;
      if (brushXselection) {
        const selection = d3.brushSelection(brushXselection.node());
        if (!range && selection) {
          /* no active selection - clear brush */
          brushXselection.call(brushX.move, null);
        } else if (range && !selection) {
          /* there is an active selection, but no brush - set the brush */
          const x0 = x(range[0]);
          const x1 = x(range[1]);
          brushXselection.call(brushX.move, [x0, x1]);
        } else if (range && selection) {
          /* there is an active selection and a brush - make sure they match */
          const moveDeltaThreshold = 1;
          const x0 = x(range[0]);
          const x1 = x(range[1]);
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

      if (d3.event.selection) {
        dispatch({
          type,
          selection: field,
          continuousNamespace: {
            isObs,
            isUserDefined,
            isDiffExp
          },
          range: [x(d3.event.selection[0]), x(d3.event.selection[1])]
        });
      } else {
        dispatch({
          type,
          selection: field,
          continuousNamespace: {
            isObs,
            isUserDefined,
            isDiffExp
          },
          range: null
        });
      }
    };
  }

  onBrushEnd(selection, x) {
    return () => {
      const { dispatch, field, isObs, isUserDefined, isDiffExp } = this.props;
      const { brushXselection } = this.state;
      const minAllowedBrushSize = 10;
      const smallAmountToAvoidInfiniteLoop = 0.1;

      // ignore programmatically generated events
      if (!d3.event.sourceEvent) return;

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

          d3.event.target.move(brushXselection, [
            d3.event.selection[0],
            procedurallyResizedBrushWidth
          ]);
        }

        dispatch({
          type: "continuous metadata histogram end",
          selection: field,
          continuousNamespace: {
            isObs,
            isUserDefined,
            isDiffExp
          },
          range: _range
        });
      } else {
        dispatch({
          type: "continuous metadata histogram cancel",
          selection: field,
          continuousNamespace: {
            isObs,
            isUserDefined,
            isDiffExp
          }
        });
      }
    };
  }

  drawHistogram(svgRef) {
    const { obsAnnotations, field, ranges } = this.props;
    const histogramCache = this.calcHistogramCache(
      obsAnnotations,
      field,
      ranges.min,
      ranges.max
    );

    const { x, y, bins, numValues } = histogramCache;

    this._histogram = { x, y, bins, numValues, svgRef };
  }

  handleColorAction() {
    const { obsAnnotations, dispatch, field, world, ranges } = this.props;

    if (obsAnnotations.hasCol(field)) {
      dispatch({
        type: "color by continuous metadata",
        colorAccessor: field,
        rangeForColorAccessor: ranges
      });
    } else if (world.varData.hasCol(field)) {
      dispatch(actions.requestSingleGeneExpressionCountsForColoringPOST(field));
    }
  }

  removeHistogram() {
    const {
      dispatch,
      field,
      colorAccessor,
      scatterplotXXaccessor,
      scatterplotYYaccessor
    } = this.props;
    dispatch({
      type: "clear user defined gene",
      data: field
    });
    if (field === colorAccessor) {
      dispatch({
        type: "reset colorscale"
      });
    }
    if (field === scatterplotXXaccessor) {
      dispatch({
        type: "set scatterplot x",
        data: null
      });
    }
    if (field === scatterplotYYaccessor) {
      dispatch({
        type: "set scatterplot y",
        data: null
      });
    }
  }

  handleSetGeneAsScatterplotX() {
    return () => {
      const { dispatch, field } = this.props;
      dispatch({
        type: "set scatterplot x",
        data: field
      });
    };
  }

  handleSetGeneAsScatterplotY() {
    return () => {
      const { dispatch, field } = this.props;
      dispatch({
        type: "set scatterplot y",
        data: field
      });
    };
  }

  renderAxesBrushBins(x, y, bins, numValues, svgRef, field) {
    /* Remove everything */
    d3.select(svgRef)
      .selectAll("*")
      .remove();

    /* BINS */
    d3.select(svgRef)
      .insert("g", "*")
      .attr("fill", "#bbb")
      .selectAll("rect")
      .data(bins)
      .enter()
      .append("rect")
      .attr("class", "bar")
      .attr("x", d => x(d.x0) + 1)
      .attr("y", d => y(d.length / numValues))
      .attr("width", d => Math.abs(x(d.x1) - x(d.x0) - 1))
      .attr("height", d => y(0) - y(d.length / numValues));

    /* BRUSH */
    const brushX = d3
      .brushX()
      /*
      emit start so that the Undoable history can save an undo point
      upon drag start, and ignore the subsequent intermediate drag events.
      */
      .on("start", this.onBrush(field, x.invert, "start").bind(this))
      .on("brush", this.onBrush(field, x.invert, "brush").bind(this))
      .on("end", this.onBrushEnd(field, x.invert).bind(this));
    const brushXselection = d3
      .select(svgRef)
      .append("g")
      .attr("class", "brush")
      .attr("data-testid", `${svgRef.dataset.testid}-brush`)
      .call(brushX);

    /* AXIS */
    d3.select(svgRef)
      .append("g")
      .attr("class", "axis axis--x")
      .attr("transform", `translate(0,${this.height - this.marginBottom})`)
      .call(d3.axisBottom(x).ticks(5));

    d3.select(svgRef)
      .selectAll(".axis--x text")
      .style("fill", "rgb(80,80,80)");

    d3.select(svgRef)
      .selectAll(".axis--x path")
      .style("stroke", "rgb(230,230,230)");

    d3.select(svgRef)
      .selectAll(".axis--x line")
      .style("stroke", "rgb(230,230,230)");

    this.setState({ brushX, brushXselection });
  }

  render() {
    const {
      field,
      colorAccessor,
      isUserDefined,
      isDiffExp,
      logFoldChange,
      pval,
      pvalAdj,
      scatterplotXXaccessor,
      scatterplotYYaccessor,
      zebra
    } = this.props;
    const field_for_id = field.replace(/\s/g, "_");
    return (
      <div
        id={`histogram_${field_for_id}`}
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
          backgroundColor: zebra ? globals.lightestGrey : "white"
        }}
      >
        <div style={{ display: "flex", justifyContent: "flex-end" }}>
          {isDiffExp || isUserDefined ? (
            <span>
              <span
                style={{ marginRight: 7 }}
                className="bp3-icon-standard bp3-icon-scatter-plot"
              />
              <ButtonGroup style={{ marginRight: 7 }}>
                <Button
                  data-testid={`plot-x-${field}`}
                  onClick={this.handleSetGeneAsScatterplotX(field).bind(this)}
                  active={scatterplotXXaccessor === field}
                  intent={scatterplotXXaccessor === field ? "primary" : "none"}
                >
                  plot x
                </Button>
                <Button
                  data-testid={`plot-y-${field}`}
                  onClick={this.handleSetGeneAsScatterplotY(field).bind(this)}
                  active={scatterplotYYaccessor === field}
                  intent={scatterplotYYaccessor === field ? "primary" : "none"}
                >
                  plot y
                </Button>
              </ButtonGroup>
            </span>
          ) : null}
          {isUserDefined ? (
            <Button
              minimal
              onClick={this.removeHistogram.bind(this)}
              style={{
                color: globals.blue,
                cursor: "pointer",
                marginLeft: 7
              }}
            >
              remove
            </Button>
          ) : null}
          <Tooltip content="Use as color scale" position="bottom">
            <Button
              onClick={this.handleColorAction.bind(this)}
              active={colorAccessor === field}
              intent={colorAccessor === field ? "primary" : "none"}
              data-testclass="colorby"
              data-testid={`colorby-${field}`}
              icon="tint"
            />
          </Tooltip>
        </div>
        <svg
          width={this.width}
          height={this.height}
          id={`histogram_${field_for_id}_svg`}
          data-testclass="histogram-plot"
          data-testid={`histogram-${field}-plot`}
          ref={svgRef => {
            this.drawHistogram(svgRef);
          }}
        />
        <div
          style={{
            display: "flex",
            justifyContent: "center"
          }}
        >
          <span
            data-testclass="brushable-histogram-field-name"
            style={{ fontStyle: "italic" }}
          >
            {field}
          </span>
        </div>

        {isDiffExp ? (
          <div
            style={{
              display: "flex",
              justifyContent: "center",
              alignItems: "baseline"
            }}
          >
            <span>
              <strong>log fold change:</strong>
              {` ${logFoldChange.toPrecision(4)}`}
            </span>
            <span
              style={{
                marginLeft: 7,
                padding: 2
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
