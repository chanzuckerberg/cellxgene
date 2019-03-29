const GraphSelection = (state = {}, action) => {
  switch (action.type) {
    case "reset World to eq Universe": {
      return {
        mode: "all"
      };
    }

    case "graph brush selection change": {
      const { brushCoords } = action;
      return {
        mode: "rect",
        brushCoords
      };
    }

    case "graph lasso selection": {
      const { polygon } = action;
      return {
        mode: "lasso",
        polygon
      };
    }

    case "graph lasso deselect":
    case "graph brush deselect": {
      return {
        mode: "all"
      };
    }

    default: {
      return state;
    }
  }
};

export default GraphSelection;
