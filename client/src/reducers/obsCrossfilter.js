/*
Reducer for the obsCrossfilter
*/

const ObsCrossfilter = (state = null, action) => {
  switch (action.type) {
    /*
		many, many actions imply a new crossfilter.  It might be simpler to
		just have a test for the existance of `obsCrossfilter` in the action,
		and assume it changes the crossfilter?
		*/
    case "annoMatrix: init complete":
    case "continuous metadata histogram start":
    case "continuous metadata histogram brush":
    case "continuous metadata histogram cancel":
    case "continuous metadata histogram end":
    case "categorical metadata filter select":
    case "categorical metadata filter deselect":
    case "categorical metadata filter none of these":
    case "categorical metadata filter all of these":
    case "graph brush change":
    case "graph brush end":
    case "graph brush cancel":
    case "graph brush deselect":
    case "graph lasso end":
    case "graph lasso cancel":
    case "graph lasso deselect":
    case "set clip quantiles":
    case "subset to selection":
    case "reset subset": {
      const { obsCrossfilter } = action;
      if (!obsCrossfilter) return state;
      return obsCrossfilter;
    }

    default: {
      return state;
    }
  }
};

export default ObsCrossfilter;
