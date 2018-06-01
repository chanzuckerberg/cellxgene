// jshint esversion: 6
const Expression = (
  state = {
    data: null,
    loading: null,
    error: null
  },
  action
) => {
  switch (action.type) {
    case "get expression started":
      return Object.assign({}, state, {
        loading: true,
        error: null
      });
    case "get expression success":
      return Object.assign({}, state, {
        error: null,
        loading: false,
        data: action.data.data
      });
    case "get expression error":
      return Object.assign({}, state, {
        data: null,
        loading: false,
        error: action.data
      });
    default:
      return state;
  }
};

export default Expression;
