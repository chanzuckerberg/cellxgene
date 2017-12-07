import { combineReducers, createStore, applyMiddleware } from 'redux';
import updateURLMiddleware from "../middleware/updateURLMiddleware";
import updateCellSelectionMiddleware from "../middleware/updateCellSelectionMiddleware";

import thunk from "redux-thunk";

import initialize from "./initialize";
import cells from "./cells";
import expression from "./expression";
import controls from "./controls";
import selectedMetadata from "./selectedMetadata";

const Reducer = combineReducers({
  initialize,
  cells,
  expression,
  controls,
  selectedMetadata,
})

let store = createStore(
  Reducer,
  window.__REDUX_DEVTOOLS_EXTENSION__ && window.__REDUX_DEVTOOLS_EXTENSION__(),
  applyMiddleware(
    thunk,
    updateURLMiddleware,
    updateCellSelectionMiddleware
  )
);

export default store;
