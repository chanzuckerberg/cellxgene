/*
Reducer for the obsCrossfilter
*/

const ObsCrossfilter = (state = null, action) => {
  if (action.obsCrossfilter) {
    return action.obsCrossfilter;
  }
  return state;
};

export default ObsCrossfilter;
