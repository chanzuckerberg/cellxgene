import _ from "lodash";

const Controls = (state = {
  _ranges: null, /* this comes from initialize, this is universe */
  allCellsOnClient: null, /* this comes from cells endpoint, this is world */
  currentCellSelection: null, /* this comes from user actions, this is current cell selection */
  graphMap: null,
  colorAccessor: null,
  colorScale: null,
  graphBrushSelection: null,
  continuousSelection: null,
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
    const graphMap = {};
    const currentCellSelection = action.data.data.metadata.slice(0);
    _.each(action.data.data.graph, (g) => { graphMap[g[0]] = [g[1], g[2]] });
    _.each(currentCellSelection, (cell) => {
      cell["__selected__"] = true;
      cell["__color__"] = "rgba(255,0,0,1)" /* initial color for all cells in all charts */
    });

    return Object.assign({}, state, {
      allCellsOnClient: action.data.data,
      currentCellSelection,
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
      continuousSelection: action.data,
      currentCellSelection: action.newSelection /* this comes from middleware */
    });
  }
  case "graph brush selection change":
    return Object.assign({}, state, {
      graphBrushSelection: action.brushCoords, /* this has already been applied in middleware but store it for next time */
      currentCellSelection: action.newSelection /* this comes from middleware */
    })
  case "graph brush deselect":
    return Object.assign({}, state, {
      graphBrushSelection: null,
      currentCellSelection: action.newSelection
    })
  case "color by continuous metadata":
    return Object.assign({}, state, {
      colorAccessor: action.colorAccessor,
      currentCellSelection: action.currentSelectionWithUpdatedColors /* this comes from middleware */
    });
  case "color by expression":
    return Object.assign({}, state, {
      colorAccessor: action.gene,
      currentCellSelection: action.currentSelectionWithUpdatedColors /* this comes from middleware */
    })
  default:
    return state;
  }
};

export default Controls;
