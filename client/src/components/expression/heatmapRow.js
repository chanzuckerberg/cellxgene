import React from "react";
import { connect } from "react-redux";
import { FaPaintBrush } from "react-icons/fa";
import HeatmapSquare from "./heatmapSquare";
import * as globals from "../../globals";
import actions from "../../actions";

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

export default HeatmapRow;
