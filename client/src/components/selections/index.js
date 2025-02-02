import React, { useRef, useEffect } from "react";
import { connect } from "react-redux";
import { H4 } from "@blueprintjs/core";
import cls from "./index.css";

import actions from "../../actions";

import { Selection } from "./selection";

const Selections = ({ selections, layoutChoice, dispatch }) => {
  const ref = useRef(null);

  useEffect(() => {
    if (ref.current) {
      ref.current.scrollTop = ref.current.scrollHeight;
    }
  }, [selections]);

  const onMouseOver = (name) => {
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
  };

  // todo: 与undo功能冲突，在确定需要该功能时再修复
  // const onMouseLeave = () => {
  //   dispatch(actions.graphLassoDeselectAction(layoutChoice.current));
  // }

  return (
    <div className={cls.root}>
      <H4>Region Sets ({layoutChoice.current})</H4>
      <div className={cls.list} ref={ref}>
        {selections.length === 0 ? (
          <span className={cls.empty}>Empty</span>
        ) : null}
        {selections.map((s) => (
          <Selection
            key={s.name}
            name={s.name}
            indexes={s.indexes}
            onMouseOver={() => onMouseOver(s.name)}
            // onMouseLeave={() => onMouseLeave()}
          />
        ))}
      </div>
    </div>
  );
};

const mapStateToProps = (state) => ({
  selections: state.graphSelection.selections,
  layoutChoice: state.layoutChoice,
});

export default connect(mapStateToProps)(Selections);
