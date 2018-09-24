// jshint esversion: 6
import React from "react";
import _ from "lodash";
import { connect } from "react-redux";
import * as globals from "../../globals";
import actions from "../../actions";
import CellSetButton from "./cellSetButtons";

@connect(state => {
  return {
    differential: state.differential,
    world: state.controls.world,
    crossfilter: _.get(state.controls, "crossfilter", null),
    selectionUpdate: _.get(state.controls, "crossfilter.updateTime", null)
  };
})
class Expression extends React.Component {
  constructor(props) {
    super(props);
    this.state = {};
  }

  handleClick(gene) {
    return () => {
      this.props.dispatch({
        type: "color by expression",
        gene: gene
      });
    };
  }

  computeDiffExp() {
    this.props.dispatch(
      actions.requestDifferentialExpression(
        this.props.differential.celllist1,
        this.props.differential.celllist2
      )
    );
  }

  render() {
    if (!this.props.differential) {
      return null;
    }
    return (
      <div style={{ marginRight: 10, marginBottom: 10 }}>
        <CellSetButton {...this.props} eitherCellSetOneOrTwo={1} />
        <CellSetButton {...this.props} eitherCellSetOneOrTwo={2} />
        <button
          style={{
            fontSize: 14,
            fontWeight: 400,
            color: "#FFF",
            padding: "0px 10px",
            height: 30,
            borderRadius: 2,
            backgroundColor:
              this.props.differential.celllist1 &&
              this.props.differential.celllist2
                ? globals.brightBlue
                : globals.lightGrey,
            border: "none",
            cursor:
              this.props.differential.celllist1 &&
              this.props.differential.celllist2
                ? "pointer"
                : "auto"
          }}
          onClick={this.computeDiffExp.bind(this)}
        >
          Compute differential expression
        </button>
      </div>
    );
  }
}

export default Expression;

// <div style={{ marginBottom: 15, width: 300 }}>
//   There are currently
//   {" " +
//     (this.props.crossfilter
//       ? this.props.crossfilter.cells.countFiltered()
//       : 0) +
//     " "}
//   cells selected, click a cell set button to store them.
// </div>
