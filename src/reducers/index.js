import { combineReducers, createStore, applyMiddleware } from 'redux';
import updateURLMiddleware from "../middleware/updateURLMiddleware";
import thunk from "redux-thunk";

import initialize from "./initialize";
import cells from "./cells";
import controls from "./controls";
import selectedMetadata from "./selectedMetadata";

const Reducer = combineReducers({
  initialize,
  cells,
  controls,
  selectedMetadata,
})

let store = createStore(
  Reducer,
  window.__REDUX_DEVTOOLS_EXTENSION__ && window.__REDUX_DEVTOOLS_EXTENSION__(),
  applyMiddleware(thunk, updateURLMiddleware)
);

export default store;
