const Dotplot = (
  state = {
    row: null,
    column: null,
  },
  action
) => {
  switch (action.type) {
    case "set dotplot row":
      console.log("in dotplot reducer row is", action.data);
      return {
        ...state,
        row: action.data,
      };
    case "set dotplot column":
      console.log("in dotplot reducer COLUMN is", action.data);

      return {
        ...state,
        column: action.data,
      };
    default:
      return state;
  }
};

export default Dotplot;
