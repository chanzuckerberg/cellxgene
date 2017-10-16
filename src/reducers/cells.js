const Cells = (state = {
  cells: null,
  loading: null,
  error: null,
}, action) => {
  switch (action.type) {
  case "request cells started":
    return Object.assign({}, state, {
      loading: true,
      error: null
    });
  case "request cells success":
    return Object.assign({}, state, {
      error: null,
      loading: false,
      cells: action.data,
    });
  case "request cells error":
    return Object.assign({}, state, {
      cells: null,
      loading: false,
      error: action.data
    });
  default:
    return state;
  }
};

export default Cells;
