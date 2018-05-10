// jshint esversion: 6
import React from "react";
import _ from "lodash";
import { connect } from "react-redux";
import * as globals from "../../globals";
import actions from "../../actions";
import CellSetButton from "./cellSetButtons";

@connect(state => {
  return {
    currentCellSelection: state.controls.currentCellSelection,
    differential: state.differential
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
      <div>
        <div style={{ margin: 10 }}>
          <div style={{ marginBottom: 15, width: 300 }}>
            There are currently
            {" " +
              _.filter(this.props.currentCellSelection, "__selected__").length +
              " "}
            cells selected, click a cell set button to store them.
          </div>
          <CellSetButton {...this.props} eitherCellSetOneOrTwo={1} />
          <CellSetButton {...this.props} eitherCellSetOneOrTwo={2} />
        </div>
        <div>
          {this.props.differential.celllist1 &&
          this.props.differential.celllist2 ? (
            <button
              style={{
                fontSize: 18,
                margin: 10,
                fontWeight: 700,
                color: "#FFF",
                padding: "12px 20px",
                backgroundColor: globals.brightBlue,
                border: "none",
                cursor: "pointer"
              }}
              onClick={this.computeDiffExp.bind(this)}
            >
              Compute differential expression
            </button>
          ) : (
            <button
              style={{
                fontSize: 18,
                margin: 10,
                fontWeight: 700,
                color: "#FFF",
                padding: "12px 20px",
                backgroundColor: globals.mediumGrey,
                border: "none"
              }}
            >
              Compute differential expression
            </button>
          )}
        </div>
      </div>
    );
  }
}

export default Expression;
