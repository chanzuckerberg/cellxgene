const Controls = (state = {
  color: null,
  continuousSelection: null,
  graphBrushSelection: null
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
  case "graph brush selection change":
    return Object.assign({}, state, {
      graphBrushSelection: action.brushCoords
    })
  case "graph brush deselect":
    return Object.assign({}, state, {
      graphBrushSelection: null
    })
  default:
    return state;
  }
};

export default Controls;
