const GraphSelection = (
  state = {
    tool: "lasso", // what selection tool mode (lasso, brush, ...)
    selection: { mode: "all" } // current selection, which is tool specific
  },
  action
) => {
  switch (action.type) {
    case "reset World to eq Universe": {
      return {
        ...state,
        selection: {
          mode: "all"
        }
      };
    }

    case "graph brush end":
    case "graph brush change": {
      const { brushCoords, sourceCoords } = action;
      return {
        ...state,
        selection: {
          mode: "within-rect",
          brushCoords,
          sourceCoords
        }
      };
    }

    case "graph lasso end": {
      const { polygon, source } = action;
      return {
        ...state,
        selection: {
          mode: "within-polygon",
          polygon,
          source
        }
      };
    }

    case "graph lasso cancel":
    case "graph brush cancel":
    case "graph lasso deselect":
    case "graph brush deselect": {
      return {
        ...state,
        selection: {
          mode: "all"
        }
      };
    }

    default: {
      return state;
    }
  }
};

export default GraphSelection;
