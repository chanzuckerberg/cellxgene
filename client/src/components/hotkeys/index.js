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
        label:
          "Set the current selection and its inverse to cell sets 1 and 2, respectively.",
        onKeyDown: async () => {
          await dispatch(actions.setCellSetFromSelection(1));
          await dispatch(actions.selectInverseSelectionAction());
          await dispatch(actions.setCellSetFromSelection(2));
          await dispatch(actions.selectInverseSelectionAction());
        },
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
export const GenesetHotkeys = (props) => {
  const { dispatch, genesets } = props;
  const inputRef = createRef();
  const hotkeys = useMemo(
    () => [
      {
        combo: "SHIFT+Q",
        global: true,
        label: "Delete the most recent geneset.",
        onKeyDown: () =>
          dispatch(
            actions.genesetDelete(Array.from(genesets.values())[0].genesetName)
          ),
      },
    ],
    [genesets]
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
