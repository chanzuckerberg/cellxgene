// jshint esversion: 6
import React from "react";
import _ from "lodash";
import { Button, Tooltip } from "@blueprintjs/core";
import { connect } from "react-redux";

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
      <Tooltip
        content="Select a group of cells, then click here to store them"
        position="top"
      >
        <Button
          style={{ marginRight: 10 }}
          type="button"
          onClick={this.set.bind(this)}
        >
          {eitherCellSetOneOrTwo}
          {": "}
          {differential[cellListName]
            ? `${differential[cellListName].length} cells`
            : "0 cells"}
        </Button>
      </Tooltip>
    );
  }
}

export default CellSetButton;
