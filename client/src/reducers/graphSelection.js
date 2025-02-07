// 选区
// 如果children不为空，则视为选区组
// interface Selection {
//   id: string; // 选区ID（不能重复）
//   emb: string; // 视图名称
//   name: string; // 选区名称 （支持自定义）
//   polygon: vec2[]; // 选区顶点坐标
//   indexes: number[]; // 选区中细胞索引
//   children: Selection[]; // 子选区
// }

import { load, save } from "../util/storage";

// import { nanoid } from 'nanoid'

// const SELECTIONS_LOCAL_STORAGE_KEY = "GRAPH_SELECTIONS"

const DEFAULT_SELECTION_GROUP_NAME = "ungrouped"

const SELECTION_LAST_ID = "selection-last-id"

const createSelection = (() => {
  const getId = (() => {
    let id = load(SELECTION_LAST_ID, -1);

    return () => {
      id += 1;
      save(SELECTION_LAST_ID, id)
      return id;
    }
  })()

  return (s) => {
    const id = `region${getId()}`;
    const name = s.name || id;

    return ({
      id,
      emb: s.emb,
      name,
      polygon: s.polygon || [],
      indexes: s.indexes || [],
      children: [],
    })
  };
})();

export const getSelectionPath = (selections, targetId, paths = []) => {
  for (let i in selections) {
    const s = selections[i]
    if (s.id === targetId) {
      return [...paths, i]
    }

    const res = getSelectionPath(s.children, targetId, [...paths, i])
    if (res.length > paths.length) {
      return res
    }
  }
  return []
}

export const getSelectionByPath = (selections, paths = []) => {
  if (paths.length === 0) {
    return selections;
  }

  const [first, ...rest] = paths;
  let target = selections[first];

  if (!target) {
    return undefined;
  }

  for (let path of rest) {
    target = target.children[path];
    if (!target) {
      return undefined;
    }
  }

  return target;
}

// 移动选区(不考虑嵌套)
// 会直接修改传入的selections
const moveSelection = (selections, fromSelectionId, toGroupId) => {
  const fromPaths = getSelectionPath(selections, fromSelectionId)
  let toPaths = getSelectionPath(selections, toGroupId)

  // 起始点必须为二级
  // 目标必须大于等于一级
  if (fromPaths.length !== 2 || toPaths.length === 0) {
    return;
  }

  // 目标路径减少到一级
  toPaths = [toPaths[0]]
  // 起始目标上级目录路径
  const fromParentPaths = fromPaths.slice(0, -1)

  const fromParentSelection = getSelectionByPath(selections, fromParentPaths);
  const fromSelection = getSelectionByPath(selections, fromPaths)
  const toSelection = getSelectionByPath(selections, toPaths)

  fromParentSelection.children = fromParentSelection.children.filter(s => s.id !== fromSelection.id)
  toSelection.children = [...toSelection.children, fromSelection]
}

// 删除选区
const removeSelection = (selections, selectionId) => {
  const paths = getSelectionPath(selections, selectionId);
  if (paths.length === 1) {
    selections.splice(paths[0], 1)
    return;
  }
  const parentPaths = paths.slice(0, -1);
  const parentSelection = getSelectionByPath(selections, parentPaths);
  parentSelection.children = parentSelection.children.filter(s => s.id !== selectionId)
}

// 修改选区
const updateSelection = (selections, selectionId, props = {}) => {
  const paths = getSelectionPath(selections, selectionId);
  const selection = getSelectionByPath(selections, paths);
  Object.entries(props).forEach(([k, v]) => {
    selection[k] = v;
  })
}

// const loadedSelections = load(SELECTIONS_LOCAL_STORAGE_KEY, [])

const GraphSelection = (
  state = {
    tool: "lasso", // what selection tool mode (lasso, brush, ...)
    selection: { mode: "all" }, // current selection, which is tool specific
    // selections: loadedSelections, // 记录所有选区
    selections: [], // 记录所有选区
  },
  action
) => {
  switch (action.type) {
    case "set clip quantiles":
    case "subset to selection":
    case "reset subset":
    case "set layout choice": {
      return {
        ...state,
        selection: {
          mode: "all",
        },
      };
    }

    case "graph brush end":
    case "graph brush change": {
      const { brushCoords } = action;
      return {
        ...state,
        selection: {
          mode: "within-rect",
          brushCoords,
        },
      };
    }

    case "graph lasso end": {
      const { emb, polygon, indexes, recordable = true } = action;

      const newState = {
        ...state,
        selection: {
          mode: "within-polygon",
          polygon,
        },
        selections: [...state.selections],
      };

      if (recordable) {
        // 新增的选区默认添加到分组(name="ungrouped")
        // Todo: 支持默认添加到其他分组，有需求再实现

        // 找到默认分组
        let group = newState.selections.find(s => s.name === DEFAULT_SELECTION_GROUP_NAME)

        if (!group) {
          // 自动创建默认分组
          group = createSelection({
            emb,
            name: DEFAULT_SELECTION_GROUP_NAME,
            polygon: [],
            indexes: [],
          })
          newState.selections.push(group)
        }

        // 创建选区
        const selection = createSelection({
          emb,
          name: undefined,
          polygon,
          indexes,
        });

        // 添加到默认分组
        group.children.push(selection);
      }

      // save(SELECTIONS_LOCAL_STORAGE_KEY, newState.selections)

      return newState;
    }

    case "graph create selection": {
      const { emb, name } = action;
      const s = createSelection({
        emb,
        name,
      })

      const newSelections = [
        ...state.selections,
        s,
      ]

      // save(SELECTIONS_LOCAL_STORAGE_KEY, newSelections)

      return ({
        ...state,
        selections: newSelections
      })
    }

    case "graph move selection": {
      const { fromId, toId } = action;

      const selections = [...state.selections]

      moveSelection(selections, fromId, toId);

      // save(SELECTIONS_LOCAL_STORAGE_KEY, selections)

      return ({
        ...state,
        selections,
      })

    }

    case "graph remove selection": {
      const { selectionId } = action;

      const selections = [...state.selections]

      removeSelection(selections, selectionId);

      // save(SELECTIONS_LOCAL_STORAGE_KEY, selections)

      return ({
        ...state,
        selections,
      })
    }

    case "graph update selection": {
      const { selectionId, props } = action;

      const selections = [...state.selections]

      updateSelection(selections, selectionId, props);

      // save(SELECTIONS_LOCAL_STORAGE_KEY, selections)

      return ({
        ...state,
        selections,
      })
    }

    case "graph lasso cancel":
    case "graph brush cancel":
    case "graph lasso deselect":
    case "graph brush deselect": {
      return {
        ...state,
        selection: {
          mode: "all",
        },
      };
    }

    default: {
      return state;
    }
  }
};

export default GraphSelection;
