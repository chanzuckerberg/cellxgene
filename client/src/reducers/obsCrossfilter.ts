/*
Reducer for the obsCrossfilter
*/

const ObsCrossfilter = (state = null, action: any) => {
  if (action.obsCrossfilter) {
    return action.obsCrossfilter;
  }
  return state;
};

export default ObsCrossfilter;
