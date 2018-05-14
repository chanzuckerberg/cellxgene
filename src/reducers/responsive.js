// jshint esversion: 6
const Responsive = (
  state = {
    width: null,
    height: null
  },
  action
) => {
  switch (action.type) {
    case "window resize":
      return Object.assign({}, state, {
        width: action.data.width,
        height: action.data.height
      });
    default:
      return state;
  }
};

export default Responsive;
