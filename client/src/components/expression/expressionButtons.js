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
  crossfilter: _.get(state.controls, "crossfilter", null),
  selectionUpdate: _.get(state.controls, "crossfilter.updateTime", null)
}))
class Expression extends React.Component {
  constructor(props) {
    super(props);
    this.state = {};
  }

  handleClick(gene) {
    const { dispatch } = this.props;
    return () => {
      dispatch({
        type: "color by expression",
        gene
      });
    };
  }

  computeDiffExp() {
    const { dispatch, differential } = this.props;
    dispatch(
      actions.requestDifferentialExpression(
        differential.celllist1,
        differential.celllist2
      )
    );
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
          Compute differential expression
        </button>
      </div>
    );
  }
}

export default Expression;
