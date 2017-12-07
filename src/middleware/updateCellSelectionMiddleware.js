import uri from "urijs";
import * as globals from "../globals";

/*
  https://medium.com/@jacobp100/you-arent-using-redux-middleware-enough-94ffe991e6
  storeInstance => functionToCallWithAnActionThatWillSendItToTheNextMiddleware => actionThatDispatchWasCalledWith => valueToUseAsTheReturnValueOfTheDispatchCall
*/

const updateCellSelectionMiddleware = (store) => {
  return (next) => {
    return (action) => {
      const s = store.getState();

      /* this is a hardcoded map of the things we need to keep an eye on and update global cell selection in response to */
      const filterJustChanged =
        action.type === "continuous selection using parallel coords brushing" ||
        action.type === "graph brush selection change" ||
        action.type === "graph brush deselect"
        ;

      if (
        !filterJustChanged ||
        !s.controls.allCellsOnClient
        /* graphMap is set at the same time as allCells, so we assume it exists */
      ) {
        return next(action); /* if the cells haven't loaded or the action wasn't a filter, bail */
      }


      // metadata has cellname, and that's all we ever need (is a key to graphMap)
      let newSelection = s.controls.allCellsOnClient.metadata.slice(0);

      /*
         in plain language...

         (a) once the cells have loaded.
         (b) each time a user changes ANY control we need to update currentCellSelection
         there are two states:

         1. control state we already know about (state.foo)
         2. control states that override states we already know about (action.foo applied instead of state.foo)

      */

      if ( /* is there a 2d graph brush selection ? */
        action.type === "graph brush selection change" ||
        s.controls.graphBrushSelection
      ) {

        const graphBrushSelection = /* it exists, so is it new or old */
          action.type === "graph brush selection change" ? action.brushCoords :
          s.controls.graphBrushSelection

        _.each(newSelection, (cell, i) => {

          if (!s.controls.graphMap[cell["CellName"]]) { return } /* this means the cells for which we do not have a graph location will be ALWAYS SELECTED */

          const coords = s.controls.graphMap[cell["CellName"]]; // [0.08005009151334168, 0.6907652173913044]

          const pointIsInsideBrushBounds = (
            globals.graphXScale(coords[0]) >= graphBrushSelection.northwestX &&
            globals.graphXScale(coords[0]) <= graphBrushSelection.southeastX &&
            globals.graphYScale(coords[1]) >= graphBrushSelection.northwestY &&
            globals.graphYScale(coords[1]) <= graphBrushSelection.southeastY
          );

          if (!pointIsInsideBrushBounds) {
            newSelection[i]["__selected__"] = false;
          } else {
            newSelection[i]["__selected__"] = true;
          }
        })
      }
      // if (
      //   action.type === "continuous selection using parallel coords brushing"
      // ) {
      //   newSelection = newSelection.filter((d) => {
      //     /* this is iterating over the enter dataset */
      //     if (action.data.every((active) => {
      //         var dim = active.dimension;
      //         // test if point is within extents for each active brush
      //         return dim.type.within(d[dim.key], active.extent, dim);
      //       })) {
      //       return true;
      //     }
      //   });
      // }

      let modifiedAction = Object.assign({}, action, {newSelection})

      return next(modifiedAction);
    }
  }
};

export default updateCellSelectionMiddleware;
