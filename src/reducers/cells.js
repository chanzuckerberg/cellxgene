const Cells = (state = {
  cells: null,
  loading: null,
  error: null,
}, action) => {
  switch (action.type) {
  case "REQUEST_CELLS":
    return Object.assign({}, state, {
      loading: true,
      error: null
    });
  case "RECEIVE_CELLS":
    return Object.assign({}, state, {
      error: null,
      cells: action.data,
    });
  case "CELLS_FETCH_ERROR":
    return Object.assign({}, state, {
      cells: null,
      error: action.data
    });
  default:
    return state;
  }
};

export default Cells;
