const Differential = (
  state = {
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
    case "reset subset":
    case "subset to selection":
      return {
        ...state,
        celllist1: null,
        celllist2: null,
      };
    default:
      return state;
  }
};

export default Differential;
