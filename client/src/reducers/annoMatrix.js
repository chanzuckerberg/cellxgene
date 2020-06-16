/*
Reducer for the annoMatrix
*/

const AnnoMatrix = (state = null, action) => {
  console.log(action);
  switch (action.type) {
    case "annoMatrix: init complete":
    case "set clip quantiles":
    case "subset to selection":
    case "reset subset": {
      const { annoMatrix } = action;
      return annoMatrix;
    }
    default: {
      return state;
    }
  }
};

export default AnnoMatrix;
