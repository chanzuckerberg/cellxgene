// jshint esversion: 6

import _ from "lodash";
import { subsetAndResetGeneLists } from "../util/stateManager/controlsHelpers";

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

    singletonHover: false,
    datasetDrawer: false,
  },
  action
) => {
  /*
  For now, log anything looking like an error to the console.
  */
  if (action.error || /error/i.test(action.type)) {
    console.error(action.error);
  }

  switch (action.type) {
    case "initial data load start": {
      return { ...state, loading: true };
    }
    case "initial data load complete": {
      /* now fully loaded */
      return {
        ...state,
        loading: false,
        error: null,
        resettingInterface: false,
      };
    }
    case "reset subset": {
      const [newUserDefinedGenes, newDiffExpGenes] = subsetAndResetGeneLists(
        state
      );
      return {
        ...state,
        resettingInterface: false,
        userDefinedGenes: newUserDefinedGenes,
        diffexpGenes: newDiffExpGenes,
      };
    }
    case "subset to selection": {
      const [newUserDefinedGenes, newDiffExpGenes] = subsetAndResetGeneLists(
        state
      );
      return {
        ...state,
        loading: false,
        error: null,
        userDefinedGenes: newUserDefinedGenes,
        diffexpGenes: newDiffExpGenes,
      };
    }
    case "request user defined gene started": {
      return {
        ...state,
        userDefinedGenesLoading: true,
      };
    }
    case "request user defined gene error": {
      return {
        ...state,
        userDefinedGenesLoading: false,
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
        userDefinedGenesLoading: false,
      };
    }
    case "request differential expression success": {
      const diffexpGenes = action.data.map((v) => v[0]);
      return {
        ...state,
        diffexpGenes,
      };
    }
    case "clear differential expression": {
      return {
        ...state,
        diffexpGenes: [],
      };
    }
    case "clear user defined gene": {
      const { userDefinedGenes } = state;
      const newUserDefinedGenes = _.filter(
        userDefinedGenes,
        (d) => d !== action.data
      );
      return {
        ...state,
        userDefinedGenes: newUserDefinedGenes,
      };
    }
    case "initial data load error": {
      return {
        ...state,
        loading: false,
        error: action.error,
      };
    }

    /*******************************
             User Events
     *******************************/
    case "change graph interaction mode":
      return {
        ...state,
        graphInteractionMode: action.data,
      };
    case "change opacity deselected cells in 2d graph background":
      return {
        ...state,
        opacityForDeselectedCells: action.data,
      };
    case "increment graph render counter": {
      const c = state.graphRenderCounter + 1;
      return {
        ...state,
        graphRenderCounter: c,
      };
    }

    /*******************************
              Scatterplot
    *******************************/
    case "set scatterplot x":
      return {
        ...state,
        scatterplotXXaccessor: action.data,
      };
    case "set scatterplot y":
      return {
        ...state,
        scatterplotYYaccessor: action.data,
      };
    case "clear scatterplot":
      return {
        ...state,
        scatterplotXXaccessor: null,
        scatterplotYYaccessor: null,
      };

    /**************************
          Dataset Drawer
     **************************/
    case "singleton hover on":
      return { ...state, singletonHover: true };

    case "singleton hover off":
      return { ...state, singletonHover: false };

    case "toggle dataset drawer":
      return { ...state, datasetDrawer: !state.datasetDrawer };

    default:
      return state;
  }
};

export default Controls;
