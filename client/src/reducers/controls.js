// jshint esversion: 6

import _ from "lodash";

const Controls = (
  state = {
    // data loading flag
    loading: true,
    error: null,

    // all of the data + selection state
    userDefinedGenes: [],
    userDefinedGenesLoading: false,
    diffexpGenes: [],

    resettingInterface: false,
    graphInteractionMode: "select",
    opacityForDeselectedCells: 0.2,
    scatterplotXXaccessor: null, // just easier to read
    scatterplotYYaccessor: null,
    graphRenderCounter: 0 /* integer as <Component key={graphRenderCounter} - a change in key forces a remount */,
    __storedStateForCelllist1__: null /* will need procedural control of brush ie., brush.extent https://bl.ocks.org/micahstubbs/3cda05ca68cba260cb81 */,
    __storedStateForCelllist2__: null
  },
  action,
  nextSharedState,
  prevSharedState
) => {
  /*
  For now, log anything looking like an error to the console.
  */
  if (action.error || /error/i.test(action.type)) {
    console.error(action.error);
  }

  switch (action.type) {
    /*****************************************************
          Initialization, World/Universe management
          and data loading.
    ******************************************************/
    case "initial data load start": {
      return { ...state, loading: true };
    }
    case "universe: column load success": {
      /* 
      we are in a loading state until the following are true:
        * universe exists (if partially)
        * embeddings (obsLayout) are loaded
        * varAnnotations index is loaded
      */
      const universeExists = !!nextSharedState.universe;
      const embeddingsExist =
        !nextSharedState.universe?.obsLayout?.isEmpty() ?? false;
      const varIndex = nextSharedState.universe.schema.annotations.var.index;
      const varAnnotationsIndexExists =
        universeExists &&
        nextSharedState.universe.varAnnotations.hasCol(varIndex);
      return {
        ...state,
        loading: !(
          universeExists &&
          embeddingsExist &&
          varAnnotationsIndexExists
        )
      };
    }
    case "initial data load complete (universe exists)": {
      /* now fully loaded */
      return {
        ...state,
        loading: false,
        error: null,
        resettingInterface: false
      };
    }
    case "reset World to eq Universe": {
      return {
        ...state,
        resettingInterface: false
      };
    }
    case "set World to current selection": {
      return {
        ...state,
        loading: false,
        error: null
      };
    }
    case "request user defined gene started": {
      return {
        ...state,
        userDefinedGenesLoading: true
      };
    }
    case "request user defined gene error": {
      return {
        ...state,
        userDefinedGenesLoading: false
      };
    }
    case "request user defined gene success": {
      const { userDefinedGenes } = state;
      const _userDefinedGenes = _.uniq(
        userDefinedGenes.concat(action.data.genes)
      );
      return {
        ...state,
        userDefinedGenes: _userDefinedGenes,
        userDefinedGenesLoading: false
      };
    }
    case "request differential expression success": {
      const { world } = prevSharedState;
      const varIndexName = world.schema.annotations.var.index;
      const _diffexpGenes = [];
      action.data.forEach(d => {
        _diffexpGenes.push(world.varAnnotations.at(d[0], varIndexName));
      });
      return {
        ...state,
        diffexpGenes: _diffexpGenes
      };
    }
    case "clear differential expression": {
      return {
        ...state,
        diffexpGenes: []
      };
    }
    case "clear user defined gene": {
      const { userDefinedGenes } = state;
      const newUserDefinedGenes = _.filter(
        userDefinedGenes,
        d => d !== action.data
      );
      return {
        ...state,
        userDefinedGenes: newUserDefinedGenes
      };
    }
    case "clear all user defined genes": {
      return {
        ...state,
        userDefinedGenes: []
      };
    }
    case "expression load error":
    case "initial data load error": {
      return {
        ...state,
        loading: false,
        error: action.error
      };
    }

    /*******************************
             User Events
     *******************************/
    case "change graph interaction mode":
      return {
        ...state,
        graphInteractionMode: action.data
      };
    case "change opacity deselected cells in 2d graph background":
      return {
        ...state,
        opacityForDeselectedCells: action.data
      };
    case "increment graph render counter": {
      const c = state.graphRenderCounter + 1;
      return {
        ...state,
        graphRenderCounter: c
      };
    }
    case "interface reset started": {
      return {
        ...state,
        resettingInterface: true
      };
    }

    /*******************************
              Scatterplot
    *******************************/
    case "set scatterplot x":
      return {
        ...state,
        scatterplotXXaccessor: action.data
      };
    case "set scatterplot y":
      return {
        ...state,
        scatterplotYYaccessor: action.data
      };
    case "clear scatterplot":
      return {
        ...state,
        scatterplotXXaccessor: null,
        scatterplotYYaccessor: null
      };

    default:
      return state;
  }
};

export default Controls;
