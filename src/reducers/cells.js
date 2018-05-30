// jshint esversion: 6
const Cells = (
  state = {
    cells: null /* world */,
    loading: null,
    error: null,

    allCells: null /* this comes from cells endpoint, this is universe */
  },
  action
) => {
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
        cells: action.data, // most recently loaded cells

        /* Universe - initialize once */
        allCells: state.allCells ? state.allCells : action.data
      });
    case "request cells error":
      return Object.assign({}, state, {
        loading: false,
        error: action.data
      });
    case "reset graph":
      return Object.assign({}, state, {
        cells: state.allCells
      });
    default:
      return state;
  }
};

export default Cells;
