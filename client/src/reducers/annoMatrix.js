/*
Reducer for the annoMatrix
*/

const AnnoMatrix = (state = null, action) => {
  console.log(action);
  switch (action.type) {
    case "annoMatrix: init complete":
    case "set clip quantiles":
    case "subset to selection":
    case "reset subset":
    case "annotation: create category":
    case "annotation: category edited":
    case "annotation: delete category": {
      const { annoMatrix } = action;
      return annoMatrix;
    }

    default: {
      // XXX debugging code
      if (action.annoMatrix)
        console.error(
          "**** action with field annoMatrix ignored - likely bug ****"
        );
      return state;
    }
  }
};

export default AnnoMatrix;
