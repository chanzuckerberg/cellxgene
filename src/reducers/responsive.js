// jshint esversion: 6
const Responsive = (
  state = {
    width: null,
    height: null
  },
  action
) => {
  switch (action.type) {
    case "resize event":
      return Object.assign({}, state, {
        width: action.data,
        height: action.data
      });
    default:
      return state;
  }
};

export default Responsive;
