// jshint esversion: 6
import React from "react";
import _ from "lodash";
import { connect } from "react-redux";
import * as globals from "../../globals";

@connect()
class CellSetButton extends React.Component {
  set() {
    const { dispatch, crossfilter, eitherCellSetOneOrTwo } = this.props;
    const set = _.map(crossfilter.allFiltered(), "name");

    dispatch({
      type: `store current cell selection as differential set ${eitherCellSetOneOrTwo}`,
      data: set
    });
  }

  render() {
    const { eitherCellSetOneOrTwo, differential } = this.props;
    const cellListName = `celllist${eitherCellSetOneOrTwo}`;
    return (
      <span style={{ marginRight: 10 }}>
        <button
          type="button"
          style={{
            color: "#FFF",
            borderRadius: 2,
            padding: "0px 10px",
            height: 30,
            backgroundColor: globals.brightBlue,
            border: "none",
            cursor: "pointer"
          }}
          onClick={this.set.bind(this)}
        >
          <span style={{ fontSize: 14, fontWeight: 700 }}>
            {" "}
            {eitherCellSetOneOrTwo}
            {": "}
          </span>
          <span
            style={{
              fontFamily: "Georgia",
              fontSize: 14,
              marginLeft: 8,
              position: "relative",
              top: 0
            }}
          >
            {differential[cellListName]
              ? `${differential[cellListName].length} cells`
              : `${0} cells`}
          </span>
        </button>
      </span>
    );
  }
}

export default CellSetButton;
