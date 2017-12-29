const Differential = (state = {
  diffExp: null,
  loading: null,
  error: null,
  celllist1: null,
  celllist2: null,
}, action) => {
  switch (action.type) {
  case "request differential expression started":
    return Object.assign({}, state, {
      loading: true,
      error: null
    });
  case "request differential expression success":
    return Object.assign({}, state, {
      error: null,
      loading: false,
      diffExp: action.data,
    });
  case "request differential expression error":
    return Object.assign({}, state, {
      loading: false,
      error: action.data
    });
  case "store current cell selection as differential set 1":
    return Object.assign({}, state, {
      celllist1: action.data
    });
  case "store current cell selection as differential set 2":
    return Object.assign({}, state, {
      celllist2: action.data
    });
  default:
    return state;
  }
};

export default Differential;
