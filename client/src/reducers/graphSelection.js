// interface Selection {
//   emb: string;
//   name: string;
//   polygon: vec2[];
//   indexes: number[];
// }

const createSelection = (() => {
  const createName = (() => {
    let index = 0;
    return () => {
      index += 1;
      return `s${index}`;
    };
  })();

  return (s) => ({
    emb: s.emb,
    name: s.name || createName(),
    polygon: s.polygon,
    indexes: s.indexes,
  });
})();

const GraphSelection = (
  state = {
    tool: "lasso", // what selection tool mode (lasso, brush, ...)
    selection: { mode: "all" }, // current selection, which is tool specific
    selections: [],
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
        const selection = createSelection({
          emb,
          name: undefined,
          polygon,
          indexes,
        });
        newState.selections.push(selection);
      }

      return newState;
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
