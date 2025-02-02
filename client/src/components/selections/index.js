import React, { useRef, useEffect, useMemo, useState } from "react";
import { connect } from "react-redux";
import { H4 } from "@blueprintjs/core";
import cls from "./index.css";

import actions from "../../actions";

import { Selection } from "./selection";
import { SelectionCreator } from './creator'
import { getSelectionByPath, getSelectionPath } from '../../reducers/graphSelection'

import { Tree } from 'antd';
import { Button } from '@blueprintjs/core'

function convertSelectionToTreeData(s, depth = 0) {
  return ({
    key: s.id,
    title: depth === 0 ? s.name : <Selection emb={s.emb} name={s.name} indexes={s.indexes} />,
    // isLeaf: depth > 0,
    children: s.children.map(c => convertSelectionToTreeData(c, depth + 1))
  });
}

const Selections = ({ selections, layoutChoice, dispatch }) => {
  const ref = useRef(null);

  const [creatorVisible, setCreatorVisible] = useState(false);

  useEffect(() => {
    if (ref.current) {
      ref.current.scrollTop = ref.current.scrollHeight;
    }
  }, [selections]);

  const reselectSelection = ([key]) => {
    if (!key) {
      return;
    }
    const paths = getSelectionPath(selections, key);
    if (paths.length !== 2) {
      // 只展示二级选区
      return;
    }
    const selection = getSelectionByPath(selections, paths)

    if (selection.emb !== layoutChoice.current) {
      return;
    }

    dispatch(
      actions.graphRedisplaySelection(
        layoutChoice.current,
        selection.polygon,
        selection.indexes
      )
    );
  };

  const onSelectionCreate = ({ name }) => {
    dispatch({
      type: "graph create selection",
      emb: layoutChoice.current,
      name,
    })
  }

  const onDrop = (e) => {
    const { dragNode, node } = e;

    console.log(dragNode, node)

    dispatch({
      type: "graph move selection",
      fromId: dragNode.key,
      toId: node.key,
    })
  };

  const treeData = useMemo(() => selections.map(s => convertSelectionToTreeData(s)), [selections])

  console.log(treeData)

  // todo: 与undo功能冲突，在确定需要该功能时再修复
  // const onMouseLeave = () => {
  //   dispatch(actions.graphLassoDeselectAction(layoutChoice.current));
  // }

  return (
    <div className={cls.root}>
      <div className={cls.header}>
        <H4>Region Sets ({layoutChoice.current})</H4>
        <Button intent="primary" onClick={() => setCreatorVisible(true)}>Create new</Button>
      </div>
      {/* <div className={cls.list} ref={ref}>
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
      </div> */}
      <Tree
        draggable
        blockNode
        showLine
        icon={false}
        showIcon={false}
        onDrop={onDrop}
        treeData={treeData}
        onSelect={reselectSelection}
      />
      <SelectionCreator visible={creatorVisible} onClose={() => setCreatorVisible(false)} onCreate={onSelectionCreate} />
    </div>
  );
};

const mapStateToProps = (state) => ({
  selections: state.graphSelection.selections,
  layoutChoice: state.layoutChoice,
});

export default connect(mapStateToProps)(Selections);
