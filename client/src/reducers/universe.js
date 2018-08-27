// jshint esversion: 6

import { Universe } from "../util/stateManager";

/*
REST API 0.1 dataframe reducer.

TODO: 0.2 reducer
*/

const UniverseReducerV01 = (state = new Universe("0.1"), action) => {
  switch (action.type) {
    case "initial data load start": {
      return Object.assign(new Universe(), state, {
        loading: true,
        error: null
      });
    }

    /* /api/v0.1/initialize response */
    case "initialize success": {
      return Object.assign(new Universe(), state).initFromInitialize(
        action.data
      );
    }

    /* /api/v0.1/cells response */
    case "request cells success": {
      return Object.assign(new Universe(), state).initFromCells(action.data);
    }

    /* /api/v0.1/expression response */
    case "expression load success": {
      return Object.assign(new Universe(), state).initFromExpression(
        action.data
      );
    }

    case "initial data load complete (universe exists)": {
      return Object.assign(new Universe(), state, {
        loading: false,
        error: null
      });
    }

    case "initial data load error": {
      return Object.assign(new Universe(), state, {
        loading: false,
        error: action.data
      });
    }

    default:
      return state;
  }
};

export default UniverseReducerV01;
