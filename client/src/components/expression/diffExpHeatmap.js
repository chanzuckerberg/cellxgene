// jshint esversion: 6
import React from "react";
import _ from "lodash";
import { connect } from "react-redux";
import * as globals from "../../globals";
import styles from "./expression.css";
import SectionHeader from "../framework/sectionHeader";
import actions from "../../actions";
import ReactAutocomplete from "react-autocomplete"; /* http://emilebres.github.io/react-virtualized-checkbox/ */
import getContrast from "font-color-contrast"; // https://www.npmjs.com/package/font-color-contrast
import FaPaintBrush from "react-icons/lib/fa/paint-brush";
import * as d3 from "d3";
import { interpolateGreys } from "d3-scale-chromatic";

class HeatmapSquare extends React.Component {
  constructor(props) {
    super(props);
    this.state = {
      value: ""
    };
  }
  render() {
    const contrastColor = getContrast(
      this.props.backgroundColor
        .substring(4, this.props.backgroundColor.length - 1)
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
          backgroundColor: this.props.backgroundColor
        }}
      >
        {this.props.text}
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
@connect(state => {
  return {
    scatterplotXXaccessor: state.controls.scatterplotXXaccessor,
    scatterplotYYaccessor: state.controls.scatterplotYYaccessor,
    colorAccessor: state.controls.colorAccessor
  };
})
class HeatmapRow extends React.Component {
  constructor(props) {
    super(props);
    this.state = {
      value: ""
    };
  }
  handleGeneColorScaleClick(gene) {
    return () => {
      this.props.dispatch(
        actions.requestSingleGeneExpressionCountsForColoringPOST(
          this.props.gene
        )
      );
    };
  }
  handleSetGeneAsScatterplotX(gene) {
    return () => {
      this.props.dispatch({
        type: "set scatterplot x",
        data: this.props.gene
      });
    };
  }
  handleSetGeneAsScatterplotY(gene) {
    return () => {
      this.props.dispatch({
        type: "set scatterplot y",
        data: this.props.gene
      });
    };
  }
  render() {
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
            {this.props.gene}
          </span>
        </div>
        <HeatmapSquare
          backgroundColor={this.props.greyColorScale(this.props.set1exp)}
          text={this.props.set1exp}
        />
        <HeatmapSquare
          backgroundColor={this.props.greyColorScale(this.props.set2exp)}
          text={this.props.set2exp}
        />
        <span
          title={this.props.aveDiff}
          style={{
            fontSize: 12,
            marginLeft: 10,
            marginRight: 10
          }}
        >
          {this.props.aveDiff.toFixed(2)}
        </span>
        <span
          onClick={this.handleSetGeneAsScatterplotX(this.props.gene).bind(this)}
          style={{
            fontSize: 16,
            color:
              this.props.scatterplotXXaccessor === this.props.gene
                ? "white"
                : globals.brightBlue,
            cursor: "pointer",
            position: "relative",
            top: 1,
            fontWeight: 700,
            marginRight: 4,
            borderRadius: 3,
            padding: "2px 3px",
            backgroundColor:
              this.props.scatterplotXXaccessor === this.props.gene
                ? globals.brightBlue
                : "inherit"
          }}
        >
          X
        </span>
        <span
          onClick={this.handleSetGeneAsScatterplotY(this.props.gene).bind(this)}
          style={{
            fontSize: 16,
            color:
              this.props.scatterplotYYaccessor === this.props.gene
                ? "white"
                : globals.brightBlue,
            cursor: "pointer",
            position: "relative",
            top: 1,
            fontWeight: 700,
            marginRight: 4,
            borderRadius: 3,
            padding: "2px 3px",
            backgroundColor:
              this.props.scatterplotYYaccessor === this.props.gene
                ? globals.brightBlue
                : "inherit"
          }}
        >
          Y
        </span>
        <span
          onClick={this.handleGeneColorScaleClick(this.props.gene).bind(this)}
          style={{
            fontSize: 16,
            cursor: "pointer",
            position: "relative",
            marginRight: 6,
            borderRadius: 3,
            padding: "0px 2px 2px 2px",
            color:
              this.props.colorAccessor === this.props.gene
                ? "white"
                : "inherit",
            backgroundColor:
              this.props.colorAccessor === this.props.gene
                ? globals.brightBlue
                : "inherit"
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

@connect(state => {
  return {
    differential: state.differential,
    allGeneNames: state.controls.allGeneNames
  };
})
class Heatmap extends React.Component {
  constructor(props) {
    super(props);
    this.state = {
      value: ""
    };
  }
  render() {
    if (!this.props.differential.diffExp)
      return <p>Select cells & compute differential to see heatmap</p>;

    const topGenesForCellSet1 = this.props.differential.diffExp.data.celllist1;
    const topGenesForCellSet2 = this.props.differential.diffExp.data.celllist2;

    const extent = d3.extent(
      _.union(
        topGenesForCellSet1.mean_expression_cellset1,
        topGenesForCellSet1.mean_expression_cellset2,
        topGenesForCellSet2.mean_expression_cellset1,
        topGenesForCellSet2.mean_expression_cellset2
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
        {topGenesForCellSet1.topgenes.map((gene, i) => {
          return (
            <HeatmapRow
              key={gene}
              gene={gene}
              greyColorScale={greyColorScale}
              aveDiff={topGenesForCellSet1.ave_diff[i]}
              set1exp={Math.floor(
                topGenesForCellSet1.mean_expression_cellset1[i]
              )}
              set2exp={Math.floor(
                topGenesForCellSet1.mean_expression_cellset2[i]
              )}
            />
          );
        })}
        {topGenesForCellSet2.topgenes.map((gene, i) => {
          return (
            <HeatmapRow
              key={gene}
              gene={gene}
              greyColorScale={greyColorScale}
              aveDiff={topGenesForCellSet2.ave_diff[i]}
              set1exp={Math.floor(
                topGenesForCellSet2.mean_expression_cellset1[i]
              )}
              set2exp={Math.floor(
                topGenesForCellSet2.mean_expression_cellset2[i]
              )}
            />
          );
        })}
      </div>
    );
  }
}

export default Heatmap;
