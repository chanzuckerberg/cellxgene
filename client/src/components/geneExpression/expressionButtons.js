// jshint esversion: 6
import React from "react";
import _ from "lodash";
import { Button, AnchorButton, Tooltip } from "@blueprintjs/core";
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
    dispatch({
      type: "clear scatterplot"
    });
  }

  render() {
    const { differential } = this.props;
    if (!differential) {
      return null;
    }
    const haveBothCellSets =
      !!differential.celllist1 && !!differential.celllist2;
    return (
      <div
        style={{
          marginRight: 10,
          marginBottom: 10,
          paddingLeft: globals.leftSidebarSectionPadding
        }}
      >
        <CellSetButton {...this.props} eitherCellSetOneOrTwo={1} />
        <CellSetButton {...this.props} eitherCellSetOneOrTwo={2} />
        {!differential.diffExp ? (
          <Tooltip
            content="Add two cells selections, see the top 15 differentially expressed genes between them"
            position="bottom"
          >
            <AnchorButton
              style={{ marginTop: 10 }}
              disabled={!haveBothCellSets}
              intent="primary"
              data-testid="diffexp_button"
              loading={differential.loading}
              fill
              type="button"
              onClick={this.computeDiffExp.bind(this)}
            >
              Compute Differential Expression
            </AnchorButton>
          </Tooltip>
        ) : null}
        {differential.diffExp ? (
          <Tooltip
            content="Remove differentially expressed gene list and clear cell selections"
            position="bottom"
          >
            <Button
              type="button"
              fill
              style={{ marginTop: 10 }}
              intent="warning"
              onClick={this.clearDifferentialExpression.bind(this)}
            >
              Clear Differential Expression
            </Button>
          </Tooltip>
        ) : null}
      </div>
    );
  }
}

export default Expression;
