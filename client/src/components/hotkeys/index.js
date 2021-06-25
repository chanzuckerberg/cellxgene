import { useHotkeys, InputGroup } from "@blueprintjs/core";
import React, { createRef, useMemo } from "react";
import actions from "../../actions";

export const GlobalHotkeys = (props) => {
  const { dispatch } = props;
  const inputRef = createRef();
  const hotkeys = useMemo(
    () => [
      {
        combo: "SHIFT+I",
        global: true,
        label: "Invert the current selection",
        onKeyDown: () => dispatch(actions.selectInverseSelectionAction()),
      },
      {
        combo: "SHIFT+1",
        global: true,
        label: "Set current selection to cell set 1.",
        onKeyDown: () => dispatch(actions.setCellSetFromSelection(1)),
      },
      {
        combo: "SHIFT+2",
        global: true,
        label: "Set current selection to cell set 2.",
        onKeyDown: () => dispatch(actions.setCellSetFromSelection(2)),
      },
      {
        combo: "SHIFT+A",
        global: true,
        label: "Select all cells",
        onKeyDown: () => dispatch(actions.selectAll()),
      },
      {
        combo: "SHIFT+M",
        global: true,
        label:
          "Sets cell groups 1 and 2 to be the current and inverse selections, respectively. (SHIFT+1+I+2+A)",
        onKeyDown: async () => {
          await dispatch(actions.setCellSetFromSelection(1));
          await dispatch(actions.selectInverseSelectionAction());
          await dispatch(actions.setCellSetFromSelection(2));
          await dispatch(actions.selectAll());
        },
      },
    ],
    []
  );
  const { handleKeyDown, handleKeyUp } = useHotkeys(hotkeys);

  return (
    <div
      role="tab"
      tabIndex={0}
      onKeyDown={handleKeyDown}
      onKeyUp={handleKeyUp}
    >
      <InputGroup inputRef={inputRef} />
    </div>
  );
};

export const DgeHotkeys = (props) => {
  const { dispatch, differential } = props;
  const inputRef = createRef();
  const hotkeys = useMemo(
    () => [
      {
        combo: "SHIFT+D",
        global: true,
        label: "Run differential expression.",
        onKeyDown: () =>
          dispatch(
            actions.requestDifferentialExpression(
              differential.celllist1,
              differential.celllist2
            )
          ),
      },
    ],
    [differential]
  );
  const { handleKeyDown, handleKeyUp } = useHotkeys(hotkeys);

  return (
    <div
      role="tab"
      tabIndex={0}
      onKeyDown={handleKeyDown}
      onKeyUp={handleKeyUp}
    >
      <InputGroup inputRef={inputRef} />
    </div>
  );
};
