// jshint esversion: 6
import uri from "urijs";
import * as globals from "../globals";

/*
  https://medium.com/@jacobp100/you-arent-using-redux-middleware-enough-94ffe991e6
  storeInstance => functionToCallWithAnActionThatWillSendItToTheNextMiddleware => actionThatDispatchWasCalledWith => valueToUseAsTheReturnValueOfTheDispatchCall
*/

/*
  What this file does:

  1. fire a filter action anywhere in the app
  2. ** this middleware checks to see the state of all the currently selected filters, including the new one
  3. ** create updated selection from a copy of all the cells presently on the client (this may be a subset of 'all')
  4. ** append that new selection to the action so that it magically appears in the reducer just because the action was fired

  This is nice because we keep a lot of filtering business logic centralized (what it means in practice to be selected)
*/

const updateCellSelectionMiddleware = store => {
  return next => {
    return action => {
      const s = store.getState();

      /* this is a hardcoded map of the things we need to keep an eye on and update global cell selection in response to */
      const filterJustChanged =
        action.type === "continuous selection using parallel coords brushing" ||
        action.type === "continuous metadata histogram brush" ||
        action.type === "graph brush selection change" ||
        action.type === "graph brush deselect" ||
        action.type === "categorical metadata filter deselect" ||
        action.type === "categorical metadata filter select" ||
        action.type === "categorical metadata filter none of these" ||
        action.type === "categorical metadata filter all of these";

      if (
        !filterJustChanged ||
        !s.controls.allCellsOnClient
        /* graph is set at the same time as allCells, so we assume it exists */
      ) {
        return next(
          action
        ); /* if the cells haven't loaded or the action wasn't a filter, bail */
      }

      /*
        - make a FRESH copy of all of the cells
        - metadata has cellname and index, and that's all we ever need to reference cell info
      */
      let newSelection = s.controls.allCellsMetadata.slice(0);
      // _.forEach(newSelection, cell => (cell.__selected__ = true));
      for (let i = 0; i < newSelection.length; i++) {
        newSelection[i].__selected__ = true;
      }

      /*
         in plain language...

         (a) once the cells have loaded.
         (b) each time a user changes ANY control we need to update allCellsMetadata
         there are two states:

         1. control state we already know about (state.foo)
         2. control states that override states we already know about (action.foo applied instead of state.foo)

      */

      if (
        (action.type ===
          "continuous selection using parallel coords brushing" &&
          s.controls.continuousSelection) ||
        s.controls.continuousSelection
      ) {
        _.each(newSelection, (cell, i) => {
          const cellExtentsAreWithinContinuousSelectionBounds = s.controls.continuousSelection.every(
            active => {
              // test if point is within extents for each active brush
              return active.dimension.type.within(
                cell[active.dimension.key],
                active.extent,
                active.dimension
              );
            }
          );

          if (!cellExtentsAreWithinContinuousSelectionBounds) {
            newSelection[i]["__selected__"] = false;
          }
        });
      }

      let modifiedAction = Object.assign({}, action, {
        newSelection
      }); /* append the result of all the filters to the action the user just triggered */

      return next(modifiedAction);
    };
  };
};

export default updateCellSelectionMiddleware;
