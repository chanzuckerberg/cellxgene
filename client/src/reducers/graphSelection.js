const GraphSelection = (
  state = {
    tool: "lasso", // what selection tool mode (lasso, brush, ...)
    selection: { mode: "all" }, // current selection, which is tool specific
    multiselect: false,
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
      const { polygon } = action;
      return {
        ...state,
        selection: {
          mode: "within-polygon",
          polygon,
        },
      };
    }
    case "graph: lasso multi-selection on": {
      return {
        ...state,
        multiselect: true,
      };
    }
    case "graph: lasso multi-selection off": {
      return {
        ...state,
        multiselect: false,
      };
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
