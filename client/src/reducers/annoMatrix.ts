/*
Reducer for the annoMatrix
*/

// eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any -- - FIXME: disabled temporarily on migrate to TS.
const AnnoMatrix = (state = null, action: any) => {
  if (action.annoMatrix) {
    return action.annoMatrix;
  }
  return state;
};

export default AnnoMatrix;
