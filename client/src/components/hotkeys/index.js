import { useHotkeys, InputGroup } from "@blueprintjs/core";
import React, { createRef, useMemo } from "react";
import actions from "../../actions";

export const GlobalHotkeys = (props) => {
  const { dispatch } = props;
  const inputRef = createRef();
  let zPressed = false;
  const hotkeys = useMemo(
    () => [
      {
        combo: "SHIFT+I",
        global: true,
        label:
          "Set the current selection and its inverse to cell sets 1 and 2, respectively.",
        onKeyDown: () =>
          dispatch(actions.setCellsFromSelectionAndInverseAction()),
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
        combo: "Z",
        global: true,
        label: "Hold to use multiple selection lassos.",
        onKeyDown: () => {
          if (!zPressed) {
            zPressed = true;
            dispatch({ type: "graph: lasso multi-selection on" });
          }
        },
        onKeyUp: () => {
          zPressed = false;
          dispatch({ type: "graph: lasso multi-selection off" });
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
export const GenesetHotkeys = (props) => {
  const { dispatch, genesets, colorAccessor } = props;
  const inputRef = createRef();
  const hotkeys = useMemo(
    () => [
      {
        combo: "SHIFT+Q",
        global: true,
        label: "Delete the most recent geneset.",
        onKeyDown: async () => {
          const geneset = Array.from(genesets.values())[0].genesetName;
          if (colorAccessor === geneset) {
            dispatch({
              type: "color by geneset mean expression",
              geneset,
            });
          }
          dispatch(actions.genesetDelete(geneset));
        },
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
