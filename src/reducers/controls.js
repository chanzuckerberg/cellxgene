const Controls = (state = {
  color: null,
  continuousSelection: null,
}, action) => {
  switch (action.type) {
  case "color changed":
    return Object.assign({}, state, {
      color: action.data,
    });
  case "continuous selection using parallel coords brushing": {
    return Object.assign({}, state, {
      continuousSelection: action.data,
    });
  }
  default:
    return state;
  }
};

export default Controls;
