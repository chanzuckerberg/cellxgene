// jshint esversion: 6
import { combineReducers, createStore, applyMiddleware } from "redux";
import thunk from "redux-thunk";
import updateURLMiddleware from "../middleware/updateURLMiddleware";
import updateCellColors from "../middleware/updateCellColors";
import { composeWithDevTools } from "redux-devtools-extension";

import config from "./config";
import differential from "./differential";
import responsive from "./responsive";
import controls from "./controls";

const Reducer = combineReducers({
  config,
  responsive,
  controls,
  differential
});

const store = createStore(
  Reducer,
  composeWithDevTools(
    applyMiddleware(thunk, updateURLMiddleware, updateCellColors)
  )
);

export default store;
