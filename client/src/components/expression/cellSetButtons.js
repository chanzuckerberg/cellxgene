// jshint esversion: 6
import React from "react";
import _ from "lodash";
import { connect } from "react-redux";
import * as globals from "../../globals";
import styles from "../framework/buttons.css";
import actions from "../../actions";

@connect()
class CellSetButton extends React.Component {
  set() {
    const { differential } = this.props;

    const set = _.map(this.props.crossfilter.allFiltered(), "name");

    if (!differential.diffExp) {
      /* diffexp needs to be cleared before we store a new set */
      this.props.dispatch({
        type:
          "store current cell selection as differential set " +
          this.props.eitherCellSetOneOrTwo,
        data: set
      });
    }
  }

  render() {
    const { differential } = this.props;
    return (
      <span style={{ marginRight: 10 }}>
        <button
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
            {this.props.eitherCellSetOneOrTwo}
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
            {differential[
              "celllist" + this.props.eitherCellSetOneOrTwo
            ]
              ? this.props.differential[
                  "celllist" + this.props.eitherCellSetOneOrTwo
                ].length + " cells"
              : 0 + " cells"}
          </span>
        </button>
      </span>
    );
  }
}

export default CellSetButton;

//
// <button style={{
//   marginLeft: 3,
//   color: "#FFF",
//   padding: "0px 10px",
//   height: 30,
//   backgroundColor: globals.brightBlue,
//   border: "none",
//   cursor: "pointer",
// }}>
//   <span style={{fontSize: 24, fontWeight: 700}}> X </span>
//
// </button>
