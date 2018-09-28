/*
https://bl.ocks.org/mbostock/4341954
https://bl.ocks.org/mbostock/34f08d5e11952a80609169b7917d4172
https://bl.ocks.org/SpaceActuary/2f004899ea1b2bd78d6f1dbb2febf771
*/
// jshint esversion: 6
import React from "react";
import _ from "lodash";
import { connect } from "react-redux";
import FaPaintBrush from "react-icons/lib/fa/paint-brush";
import * as d3 from "d3";
import memoize from "memoize-one";
import * as globals from "../../globals";

@connect(state => ({
  initializeRanges: _.get(state.controls.world, "summary.obs"),
  colorAccessor: state.controls.colorAccessor,
  colorScale: state.controls.colorScale,
  obsAnnotations: _.get(state.controls.world, "obsAnnotations", null)
}))
class HistogramBrush extends React.Component {
  calcHistogramCache = memoize((obsAnnotations, metadataField, ranges) => {
    // recalculate expensive stuff
    const allValuesForContinuousFieldAsArray = _.map(
      obsAnnotations,
      metadataField
    );
    const histogramCache = {};

    histogramCache.x = d3
      .scaleLinear()
      .domain([ranges.min, ranges.max])
      .range([0, this.width]);

    histogramCache.y = d3
      .scaleLinear()
      .range([this.height - this.marginBottom, 0]);
    // .range([height - margin.bottom, margin.top]);

    histogramCache.bins = d3
      .histogram()
      .domain(histogramCache.x.domain())
      .thresholds(40)(allValuesForContinuousFieldAsArray);

    histogramCache.numValues = allValuesForContinuousFieldAsArray.length;

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
      const { dispatch, metadataField } = this.props;
      if (d3.event.selection) {
        dispatch({
          type: "continuous metadata histogram brush",
          selection: metadataField,
          range: [x(d3.event.selection[0]), x(d3.event.selection[1])]
        });
      } else {
        dispatch({
          type: "continuous metadata histogram brush",
          selection: metadataField,
          range: null
        });
      }
    };
  }

  drawHistogram(svgRef) {
    const { obsAnnotations, metadataField, ranges } = this.props;
    const histogramCache = this.calcHistogramCache(
      obsAnnotations,
      metadataField,
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
        .call(
          d3
            .brushX()
            .on("end", this.onBrush(metadataField, x.invert).bind(this))
        );

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
        .attr("font-weight", "bold")
        .text(metadataField);

      this.setState({ brush, xAxis });
    }
  }

  handleColorAction() {
    const { dispatch, metadataField, initializeRanges } = this.props;
    dispatch({
      type: "color by continuous metadata",
      colorAccessor: metadataField,
      rangeMaxForColorAccessor: initializeRanges[metadataField].range.max
    });
  }

  render() {
    const { metadataField, ranges, colorAccessor } = this.props;
    return (
      <div
        style={{
          marginTop: 10,
          position: "flex"
        }}
        id={`histogram_${metadataField}`}
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
            fontSize: 16,
            marginLeft: 4,
            // padding: this.props.colorAccessor === this.props.metadataField ? 3 : "auto",
            borderRadius: 3,
            color:
              colorAccessor === metadataField ? globals.brightBlue : "black",
            // backgroundColor: this.props.colorAccessor === this.props.metadataField ? globals.brightBlue : "inherit",
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
