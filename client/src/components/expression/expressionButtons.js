// jshint esversion: 6
import React from "react";
import _ from "lodash";
import { connect } from "react-redux";
import * as globals from "../../globals";
import actions from "../../actions";
import CellSetButton from "./cellSetButtons";

@connect(state => ({
  differential: state.differential,
  world: state.controls.world,
  crossfilter: state.controls.crossfilter,
  selectionUpdate: _.get(state.controls, "crossfilter.updateTime", null)
}))
class Expression extends React.Component {
  constructor(props) {
    super(props);
    this.state = {};
  }

  computeDiffExp() {
    const { dispatch, differential } = this.props;
    if (differential.celllist1 && differential.celllist2) {
      dispatch(
        actions.requestDifferentialExpression(
          differential.celllist1,
          differential.celllist2
        )
      );
    }
  }

  clearDifferentialExpression() {
    const { dispatch, differential } = this.props;
    dispatch({
      type: "clear differential expression",
      diffExp: differential.diffExp
    });
  }

  render() {
    const { differential } = this.props;
    if (!differential) {
      return null;
    }
    return (
      <div style={{ marginRight: 10, marginBottom: 10 }}>
        <CellSetButton {...this.props} eitherCellSetOneOrTwo={1} />
        <CellSetButton {...this.props} eitherCellSetOneOrTwo={2} />
        {!differential.diffExp ? (
          <button
            type="button"
            style={{
              fontSize: 14,
              fontWeight: 400,
              color: "#FFF",
              padding: "0px 10px",
              height: 30,
              borderRadius: 2,
              backgroundColor:
                differential.celllist1 && differential.celllist2
                  ? globals.brightBlue
                  : globals.lightGrey,
              border: "none",
              cursor:
                differential.celllist1 && differential.celllist2
                  ? "pointer"
                  : "auto"
            }}
            onClick={this.computeDiffExp.bind(this)}
          >
            {" "}
            Compute Differential Expression{" "}
          </button>
        ) : null}
        {differential.diffExp ? (
          <button
            type="button"
            style={{
              fontSize: 14,
              fontWeight: 400,
              color: "#FFF",
              padding: "0px 10px",
              height: 30,
              borderRadius: 2,
              backgroundColor: globals.brightBlue,
              border: "none",
              cursor:
                differential.celllist1 && differential.celllist2
                  ? "pointer"
                  : "auto"
            }}
            onClick={this.clearDifferentialExpression.bind(this)}
          >
            {" "}
            Clear Differential Expression{" "}
          </button>
        ) : null}
      </div>
    );
  }
}

export default Expression;
