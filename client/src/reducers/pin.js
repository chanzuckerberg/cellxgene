/*
Reducer for the obsCrossfilter
*/

const initialState = {
  left: false,
  right: false,
};

const PinReducer = (state = initialState, action) => {
  switch (action.type) {
    case "pin: update": {
      const { loc, pinned } = action;
      return {
        ...state,
        [loc]: pinned,
      };
    }
    default:
      return state;
  }
};

export default PinReducer;
