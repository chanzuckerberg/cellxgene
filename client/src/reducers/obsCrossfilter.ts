/*
Reducer for the obsCrossfilter
*/

// eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any -- - FIXME: disabled temporarily on migrate to TS.
const ObsCrossfilter = (state = null, action: any) => {
  if (action.obsCrossfilter) {
    return action.obsCrossfilter;
  }
  return state;
};

export default ObsCrossfilter;
