const Dotplot = (
  state = {
    row: null,
    column: null,
  },
  action
) => {
  switch (action.type) {
    case "set dotplot row":
      return {
        ...state,
        row: action.data,
      };
    case "set dotplot column":
      return {
        ...state,
        column: action.data,
      };
    case "toggle dotplot":
      return {
        ...state,
        row: null,
        column: null,
      };
    default:
      return state;
  }
};

export default Dotplot;
