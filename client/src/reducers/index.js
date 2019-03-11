// jshint esversion: 6
import { combineReducers, createStore, applyMiddleware } from "redux";
import thunk from "redux-thunk";
import { composeWithDevTools } from "redux-devtools-extension";
import updateCellColors from "../middleware/updateCellColors";

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
  composeWithDevTools(applyMiddleware(thunk, updateCellColors))
);

export default store;
