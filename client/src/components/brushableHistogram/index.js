/*
https://bl.ocks.org/mbostock/4341954
https://bl.ocks.org/mbostock/34f08d5e11952a80609169b7917d4172
https://bl.ocks.org/SpaceActuary/2f004899ea1b2bd78d6f1dbb2febf771
*/
// jshint esversion: 6
import React from "react";
import _ from "lodash";
import { connect } from "react-redux";
import { kvCache } from "../../util/stateManager/"
import FaPaintBrush from "react-icons/lib/fa/paint-brush";
import * as d3 from "d3";
import memoize from "memoize-one";
import * as globals from "../../globals";
import actions from "../../actions";

@connect(state => ({
  world: state.controls.world,
  crossfilter: state.controls.crossfilter,
  initializeRanges: _.get(state.controls.world, "summary.obs"),
  colorAccessor: state.controls.colorAccessor,
  colorScale: state.controls.colorScale,
  obsAnnotations: _.get(state.controls.world, "obsAnnotations", null)
}))
class HistogramBrush extends React.Component {
  calcHistogramCache = memoize((obsAnnotations, field, ranges) => {
    const { world } = this.props
    const histogramCache = {};

    histogramCache.y = d3
      .scaleLinear()
      .range([this.height - this.marginBottom, 0]);
    // .range([height - margin.bottom, margin.top]);

    if (obsAnnotations[0][field]) {

      // recalculate expensive stuff
      const allValuesForContinuousFieldAsArray = _.map(obsAnnotations, field);

      histogramCache.x = d3
        .scaleLinear()
        .domain([ranges.min, ranges.max])
        .range([0, this.width]);

      histogramCache.bins = d3
        .histogram()
        .domain(histogramCache.x.domain())
        .thresholds(40)(allValuesForContinuousFieldAsArray);

      histogramCache.numValues = allValuesForContinuousFieldAsArray.length;

    } else if (kvCache.get(world.varDataCache, field)) {
      /* it's not in observations, so it's a gene, but let's check to make sure */
      const varValues = kvCache.get(world.varDataCache, field)

      histogramCache.x = d3
        .scaleLinear()
        .domain(d3.extent(varValues)) /* replace this if we have ranges for genes back from server like we do for annotations on cells */
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

    this.width = 300;
    this.height = 100;
    this.marginBottom = 20;

    this.state = {
      brush: null
    };
  }

  onBrush(selection, x) {
    return () => {
      const { dispatch, field } = this.props;
      if (d3.event.selection) {
        dispatch({
          type: "continuous metadata histogram brush",
          selection: field,
          range: [x(d3.event.selection[0]), x(d3.event.selection[1])]
        });
      } else {
        dispatch({
          type: "continuous metadata histogram brush",
          selection: field,
          range: null
        });
      }
    };
  }

  drawHistogram(svgRef) {
    const { obsAnnotations, field, ranges } = this.props;
    const histogramCache = this.calcHistogramCache(
      obsAnnotations,
      field,
      ranges
    );

    const { x, y, bins, numValues } = histogramCache;
    d3.select(svgRef)
      .selectAll(".bar")
      .remove();

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

    if (!this.state.brush && !this.state.axis) {
      const brush = d3
        .select(svgRef)
        .append("g")
        .attr("class", "brush")
        .call(d3.brushX().on("end", this.onBrush(field, x.invert).bind(this)));

      const xAxis = d3
        .select(svgRef)
        .append("g")
        .attr("class", "axis axis--x")
        .attr("transform", `translate(0,${this.height - this.marginBottom})`)
        .call(d3.axisBottom(x).ticks(5))
        .append("text")
        .attr("x", this.width - 2)
        .attr("y", -6)
        .attr("fill", "#000")
        .attr("text-anchor", "end")
        .text(field);

      d3.select(svgRef)
        .selectAll(".axis--x text")
        .style("fill", "rgb(80,80,80)")
        .style("font-size", 12);

      d3.select(svgRef)
        .selectAll(".axis--x path")
        .style("stroke", "rgb(230,230,230)");

      d3.select(svgRef)
        .selectAll(".axis--x line")
        .style("stroke", "rgb(230,230,230)");

      this.setState({ brush, xAxis });
    }
  }

  handleColorAction() {
    const { obsAnnotations, ranges, dispatch, field, world, initializeRanges } = this.props;

    if (obsAnnotations[0][field]) {
      dispatch({
        type: "color by continuous metadata",
        colorAccessor: field,
        rangeMaxForColorAccessor: initializeRanges[field].range.max
      });
    } else if (kvCache.get(world.varDataCache, field)) {
      dispatch(
        actions.requestSingleGeneExpressionCountsForColoringPOST(
          field
        )
      );
    }

  }

  render() {

    const { field, ranges, colorAccessor } = this.props;

    return (
      <div
        style={{
          marginTop: 10,
          position: "flex"
        }}
        id={`histogram_${field}`}
      >
        <svg
          width={this.width}
          height={this.height}
          ref={svgRef => {
            this.drawHistogram(svgRef);
          }}
        >
          {ranges.min}
          {" to "}
          {ranges.max}
        </svg>
        <span
          onClick={this.handleColorAction.bind(this)}
          style={{
            fontSize: 14,
            marginLeft: 4,
            // padding: this.props.colorAccessor === this.props.field ? 3 : "auto",
            borderRadius: 3,
            color: colorAccessor === field ? globals.brightBlue : "black",
            // backgroundColor: this.props.colorAccessor === this.props.field ? globals.brightBlue : "inherit",
            position: "relative",
            top: -29,
            cursor: "pointer"
          }}
        >
          <FaPaintBrush />
        </span>
      </div>
    );
  }
}

export default HistogramBrush;
