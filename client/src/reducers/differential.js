// jshint esversion: 6

const Differential = (
  state = {
    diffExp: null,
    loading: null,
    error: null,
    celllist1: null,
    celllist2: null,
  },
  action
) => {
  switch (action.type) {
    case "request differential expression started":
      return {
        ...state,
        loading: true,
        error: null,
      };
    case "request differential expression success":
      return {
        ...state,
        error: null,
        loading: false,
        diffExp: action.data,
      };
    case "request differential expression error":
      return {
        ...state,
        loading: false,
        error: action.data,
      };
    case "store current cell selection as differential set 1":
      return {
        ...state,
        celllist1: action.data,
      };
    case "store current cell selection as differential set 2":
      return {
        ...state,
        celllist2: action.data,
      };
    case "clear differential expression":
      return {
        ...state,
        diffExp: null,
        celllist1: null,
        celllist2: null,
      };
    default:
      return state;
  }
};

export default Differential;
