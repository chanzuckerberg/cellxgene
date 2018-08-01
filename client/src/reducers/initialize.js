// jshint esversion: 6
const Initialize = (
  state = {
    data: null,
    loading: null,
    error: null
  },
  action
) => {
  switch (action.type) {
    case "initialize started":
      return Object.assign({}, state, {
        loading: true,
        error: null
      });
    case "initialize success":
      return Object.assign({}, state, {
        error: null,
        loading: false,
        data: action.data
      });
    case "initialize error":
      return Object.assign({}, state, {
        data: null,
        loading: false,
        error: action.data
      });
    default:
      return state;
  }
};

export default Initialize;
