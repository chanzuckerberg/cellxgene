// jshint esversion: 6
import React from "react";
import _ from "lodash";
import { connect } from "react-redux";
import * as globals from "../../globals";

@connect()
class CellSetButton extends React.Component {
  set() {
    const {
      differential,
      crossfilter,
      dispatch,
      eitherCellSetOneOrTwo
    } = this.props;

    const set = _.map(crossfilter.allFiltered(), "name");

    if (!differential.diffExp) {
      /* diffexp needs to be cleared before we store a new set */
      dispatch({
        type: `store current cell selection as differential set ${eitherCellSetOneOrTwo}`,
        data: set
      });
    }
  }

  render() {
    const { differential, eitherCellSetOneOrTwo } = this.props;
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
            backgroundColor: !differential.diffExp
              ? globals.brightBlue
              : globals.lightGrey,
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
              fontSize: 14,
              marginLeft: 8,
              position: "relative",
              top: 0
            }}
          >
            {differential[cellListName]
              ? `${differential[cellListName].length} cells`
              : "0 cells"}
          </span>
        </button>
      </span>
    );
  }
}

export default CellSetButton;
