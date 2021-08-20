/*
Reducer for the annoMatrix
*/

const AnnoMatrix = (state = null, action) => {
  if (action.annoMatrix) {
    return action.annoMatrix;
  }
  return state;
};

export default AnnoMatrix;
