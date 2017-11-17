const Controls = (state = {
  colorAccessor: null,
  colorScale: null,
  continuousSelection: null,
  graphBrushSelection: null,
  axesHaveBeenDrawn: false,
}, action) => {
  switch (action.type) {
  case "color changed":
    return Object.assign({}, state, {
      colorAccessor: action.colorAccessor,
      colorScale: d3.scaleLinear()
        .domain([0, action.rangeMaxForColorAccessor])
        .range([1,0])
    });
  case "continuous selection using parallel coords brushing": {
    return Object.assign({}, state, {
      continuousSelection: action.data,
    });
  }
  /* on load, set the selection to 'all', if reactive is true */
  case "initialize success":
    let allCellNames = null;
    if (action.data.data.reactive) { /* we have metadata, get all cell names */
      allCellNames = action.data.data.metadata.map((cell) => {
        return cell.CellName
      })
    }
    return Object.assign({}, state, {
      continuousSelection: allCellNames
    });
  case "parallel coordinates axes have been drawn":
    return Object.assign({}, state, {
      axesHaveBeenDrawn: true
    });
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
