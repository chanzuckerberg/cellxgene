import React from "react";
import { connect } from "react-redux";
import { H4 } from "@blueprintjs/core";

import * as globals from "../../globals";

import actions from "../../actions";

import { Selection } from "./selection";

@connect((state) => ({
  selections: state.graphSelection.selections,
  layoutChoice: state.layoutChoice,
}))
class Selections extends React.Component {
  constructor(props) {
    super(props);

    this.onMouseOver = this.onMouseOver.bind(this);
    // this.onMouseLeave = this.onMouseLeave.bind(this);
    this.handleSelectionsRef = this.handleSelectionsRef.bind(this);

    this.selectionsRef = null;
  }

  componentDidUpdate(prevProps) {
    const { selections: prevSelections } = prevProps;
    const { selections } = this.props;
    if (selections.length !== prevSelections.length) {
      if (this.selectionsRef) {
        this.selectionsRef.scrollTop = this.selectionsRef.scrollHeight;
      }
    }
  }

  handleSelectionsRef(ref) {
    this.selectionsRef = ref;
  }

  onMouseOver(name) {
    const { dispatch, selections, layoutChoice } = this.props;
    const sel = selections.find((s) => s.name === name);
    if (!sel) {
      return;
    }
    dispatch(
      actions.graphRedisplaySelection(
        layoutChoice.current,
        sel.polygon,
        sel.indexes
      )
    );
  }

  // todo: 与undo功能冲突，在确定需要该功能时再修复
  // onMouseLeave() {
  //   const { dispatch, layoutChoice } = this.props;
  //   dispatch(actions.graphLassoDeselectAction(layoutChoice.current));
  // }

  render() {
    const { selections, layoutChoice } = this.props;

    const currentEmbSelections = selections.filter(
      (s) => s.emb === layoutChoice.current
    );

    return (
      <div
        style={{
          padding: globals.leftSidebarSectionPadding,
          display: "flex",
          flexDirection: "column",
          alignItems: "stretch",
          rowGap: "1em",
        }}
      >
        <H4>Selections ({layoutChoice.current})</H4>
        <div
          ref={this.handleSelectionsRef}
          style={{
            display: "flex",
            flexDirection: "column",
            alignItems: "stretch",
            rowGap: "1em",
            overflowY: "auto",
            maxHeight: 200,
          }}
        >
          {currentEmbSelections.length === 0 ? (
            <span style={{ opacity: 0.5 }}>Empty</span>
          ) : null}
          {currentEmbSelections.map((s) => (
            <Selection
              key={s.name}
              name={s.name}
              indexes={s.indexes}
              onMouseOver={() => this.onMouseOver(s.name)}
              // onMouseLeave={() => this.onMouseLeave()}
            />
          ))}
        </div>
      </div>
    );
  }
}

export default Selections;
