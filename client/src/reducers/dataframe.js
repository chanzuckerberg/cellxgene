"use strict";
// jshint esversion: 6

import _ from "lodash";
import { Universe } from "../util/stateManager";

/*
REST API 0.1 dataframe reducer.

TODO: 0.2 reducer
*/

const DataFrameV01 = (state = new Universe("0.1"), action) => {
  switch (action.type) {
    case "dataframe load start": {
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

    case "dataframe load complete (universe exists)": {
      return Object.assign(new Universe(), state, {
        loading: false,
        error: null
      });
    }

    case "dataframe load error": {
      return Object.assign(new Universe(), state, {
        loading: false,
        error: action.data
      });
    }

    default:
      return state;
  }
};

export default DataFrameV01;
