import _ from "lodash";

const Controls = (state = {
  _ranges: null, /* this comes from initialize, this is universe */
  allCellsOnClient: null, /* this comes from cells endpoint, this is world */
  currentCellSelection: null, /* this comes from user actions, this is current cell selection */
  graphMap: null,
  colorAccessor: null,
  colorScale: null,
  graphBrushSelection: null,
  axesHaveBeenDrawn: false,
}, action) => {
  switch (action.type) {
  /* * * * * * * * * * * * * * * * * *
        Keep a copy of data
  * * * * * * * * * * * * * * * * * */
  case "initialize success":
    return Object.assign({}, state, {
      _ranges: action.data.data.ranges
    });
  case "request cells success":
    const graphMap = {}
    _.each(action.data.data.graph, (g) => { graphMap[g[0]] = [g[1], g[2]] })

    return Object.assign({}, state, {
      allCellsOnClient: action.data.data,
      currentCellSelection: action.data.data.metadata,
      graphMap,
      graphBrushSelection: null, /* if we are getting new cells from the server, the layout (probably? definitely?) just changed, so this is now irrelevant, and we WILL need to call a function to reset state of this kind when cells success happens */
    });
  /* * * * * * * * * * * * * * * * * *
            User events
  * * * * * * * * * * * * * * * * * */
  case "parallel coordinates axes have been drawn":
    return Object.assign({}, state, {
      axesHaveBeenDrawn: true
    });
  case "continuous selection using parallel coords brushing": {
    return Object.assign({}, state, {
      currentCellSelection: action.data,
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
  case "color changed":
    return Object.assign({}, state, {
      colorAccessor: action.colorAccessor,
      colorScale: d3.scaleLinear()
        .domain([0, action.rangeMaxForColorAccessor])
        .range([1,0])
    });
  default:
    return state;
  }
};

export default Controls;
