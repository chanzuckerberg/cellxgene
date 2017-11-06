const Controls = (state = {
  color: null,
}, action) => {
  switch (action.type) {
  case "color changed":
    return Object.assign({}, state, {
      color: action.data,
    });
  default:
    return state;
  }
};

export default Controls;
