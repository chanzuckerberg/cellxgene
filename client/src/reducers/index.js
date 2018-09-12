// jshint esversion: 6
import { combineReducers, createStore, applyMiddleware } from "redux";
import thunk from "redux-thunk";
import updateURLMiddleware from "../middleware/updateURLMiddleware";
import updateCellColors from "../middleware/updateCellColors";

import differential from "./differential";
import responsive from "./responsive";
import controls2 from "./controls2";

const Reducer = combineReducers({
  responsive,
  controls2,
  differential
});

const store = createStore(
  Reducer,
  window.__REDUX_DEVTOOLS_EXTENSION__ && window.__REDUX_DEVTOOLS_EXTENSION__(),
  applyMiddleware(thunk, updateURLMiddleware, updateCellColors)
);

export default store;
