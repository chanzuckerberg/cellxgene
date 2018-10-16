// jshint esversion: 6
import React from "react";
import _ from "lodash";
import memoize from "memoize-one";
import { connect } from "react-redux";
import ReactAutocomplete from "react-autocomplete"; /* http://emilebres.github.io/react-virtualized-checkbox/ */
import getContrast from "font-color-contrast"; // https://www.npmjs.com/package/font-color-contrast
import { FaPaintBrush } from "react-icons/fa";
import * as d3 from "d3";
import { interpolateGreys } from "d3-scale-chromatic";
import * as globals from "../../globals";
import actions from "../../actions";

class HeatmapSquare extends React.Component {
  constructor(props) {
    super(props);
    this.state = {
      value: ""
    };
  }

  render() {
    const { backgroundColor, text } = this.props;
    const contrastColor = getContrast(
      backgroundColor
        .substring(4, backgroundColor.length - 1)
        .replace(/ /g, "")
        .split(",")
    );
    return (
      <p
        style={{
          padding: "12px 6px",
          textAlign: "center",
          color: contrastColor,
          width: 40,
          flexShrink: 0,
          fontSize: 12,
          margin: 0,
          backgroundColor
        }}
      >
        {text}
      </p>
    );
  }
}

/**********************************
***********************************
***********************************
              Row
***********************************
***********************************
**********************************/
@connect(state => ({
  scatterplotXXaccessor: state.controls.scatterplotXXaccessor,
  scatterplotYYaccessor: state.controls.scatterplotYYaccessor,
  colorAccessor: state.controls.colorAccessor
}))
class HeatmapRow extends React.Component {
  constructor(props) {
    super(props);
    this.state = {
      value: ""
    };
  }

  handleGeneColorScaleClick() {
    return () => {
      const { dispatch, gene } = this.props;
      dispatch(actions.requestSingleGeneExpressionCountsForColoringPOST(gene));
    };
  }

  handleSetGeneAsScatterplotX() {
    return () => {
      const { dispatch, gene } = this.props;
      dispatch({
        type: "set scatterplot x",
        data: gene
      });
    };
  }

  handleSetGeneAsScatterplotY() {
    return () => {
      const { dispatch, gene } = this.props;
      dispatch({
        type: "set scatterplot y",
        data: gene
      });
    };
  }

  render() {
    const {
      gene,
      aveDiff,
      set1exp,
      set2exp,
      greyColorScale,
      scatterplotXXaccessor,
      scatterplotYYaccessor,
      colorAccessor
    } = this.props;
    return (
      <div
        style={{
          width: 220,
          display: "flex",
          justifyContent: "flex-start",
          alignItems: "baseline"
        }}
      >
        <div style={{ width: 150, flexShrink: 0 }}>
          <span
            style={{
              fontSize: 14
            }}
          >
            {gene}
          </span>
        </div>
        <HeatmapSquare
          backgroundColor={greyColorScale(set1exp)}
          text={set1exp.toPrecision(2)}
        />
        <HeatmapSquare
          backgroundColor={greyColorScale(set2exp)}
          text={set2exp.toPrecision(2)}
        />
        <span
          title={aveDiff}
          style={{
            fontSize: 12,
            marginLeft: 10,
            marginRight: 10
          }}
        >
          {aveDiff.toFixed(2)}
        </span>
        <span
          onClick={this.handleSetGeneAsScatterplotX(gene).bind(this)}
          style={{
            fontSize: 16,
            color:
              scatterplotXXaccessor === gene ? "white" : globals.brightBlue,
            cursor: "pointer",
            position: "relative",
            top: 1,
            fontWeight: 700,
            marginRight: 4,
            borderRadius: 3,
            padding: "2px 3px",
            backgroundColor:
              scatterplotXXaccessor === gene ? globals.brightBlue : "inherit"
          }}
        >
          X
        </span>
        <span
          onClick={this.handleSetGeneAsScatterplotY(gene).bind(this)}
          style={{
            fontSize: 16,
            color:
              scatterplotYYaccessor === gene ? "white" : globals.brightBlue,
            cursor: "pointer",
            position: "relative",
            top: 1,
            fontWeight: 700,
            marginRight: 4,
            borderRadius: 3,
            padding: "2px 3px",
            backgroundColor:
              scatterplotYYaccessor === gene ? globals.brightBlue : "inherit"
          }}
        >
          Y
        </span>
        <span
          onClick={this.handleGeneColorScaleClick(gene).bind(this)}
          style={{
            fontSize: 16,
            cursor: "pointer",
            position: "relative",
            marginRight: 6,
            borderRadius: 3,
            padding: "0px 2px 2px 2px",
            color: colorAccessor === gene ? "white" : "inherit",
            backgroundColor:
              colorAccessor === gene ? globals.brightBlue : "inherit"
          }}
        >
          <FaPaintBrush style={{ display: "inline-block" }} />
        </span>
      </div>
    );
  }
}

/**********************************
***********************************
***********************************
        DiffExp Heatmap
***********************************
***********************************
**********************************/

@connect(state => ({
  differential: state.differential,
  world: state.controls.world
}))
class Heatmap extends React.Component {
  constructor(props) {
    super(props);
    this.state = {
      value: ""
    };
  }

  // XXX TODO unused at the moment
  // getAllGeneNames = memoize(world =>
  //   _.map(this.props.world.varAnnotations, "name")
  // );

  render() {
    const { world, differential } = this.props;
    if (!differential.diffExp) {
      return <p>Select cells & compute differential to see heatmap</p>;
    }

    // summarize the information for display.
    const topGenes = _.map(differential.diffExp, val => ({
      varIndex: val[0],
      geneName: world.varAnnotations[val[0]].name,
      avgDiff: val[1],
      set1AvgExp: val[4],
      set2AvgExp: val[5]
    }));

    // average expression extent
    const extent = d3.extent(
      _.concat(
        _.map(differential.diffExp, val => val[4]),
        _.map(differential.diffExp, val => val[5])
      )
    );

    const greyColorScale = d3
      .scaleSequential()
      .domain(extent)
      .interpolator(interpolateGreys);

    return (
      <div>
        <div
          style={{
            display: "flex",
            justifyContent: "flex-start",
            width: 400,
            fontWeight: 700
          }}
        >
          <p style={{ marginRight: 110 }}>Gene</p>
          <p style={{ marginRight: 25 }}>1</p>
          <p style={{ marginRight: 20 }}>2</p>
          <p>ave diff</p>
        </div>
        {topGenes.map(g => {
          const { geneName, avgDiff, set1AvgExp, set2AvgExp } = g;
          return (
            <HeatmapRow
              key={geneName}
              gene={geneName}
              greyColorScale={greyColorScale}
              aveDiff={avgDiff}
              set1exp={set1AvgExp}
              set2exp={set2AvgExp}
            />
          );
        })}
      </div>
    );
  }
}

export default Heatmap;
