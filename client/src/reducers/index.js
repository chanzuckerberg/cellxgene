// jshint esversion: 6
import { combineReducers, createStore, applyMiddleware } from "redux";
import updateURLMiddleware from "../middleware/updateURLMiddleware";
import updateCellColors from "../middleware/updateCellColors";

import thunk from "redux-thunk";

import expression from "./expression";
import differential from "./differential";
import responsive from "./responsive";
import dataframe from "./dataframe";
import controls2 from "./controls2";

const Reducer = combineReducers({
  // XXX: old - needs refactoring
  expression,
  differential,

  // New Redux refactor
  dataframe,
  responsive,
  controls2
});

let store = createStore(
  Reducer,
  window.__REDUX_DEVTOOLS_EXTENSION__ && window.__REDUX_DEVTOOLS_EXTENSION__(),
  applyMiddleware(thunk, updateURLMiddleware, updateCellColors)
);

export default store;
