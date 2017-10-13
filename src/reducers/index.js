import { combineReducers, createStore, applyMiddleware } from 'redux';
import updateURLMiddleware from "../middleware/updateURLMiddleware";
import thunk from "redux-thunk";
import cells from "./cells";
import url from "./url";

const Reducer = combineReducers({
  cells,
  url,
})

let store = createStore(
  Reducer,
  window.__REDUX_DEVTOOLS_EXTENSION__ && window.__REDUX_DEVTOOLS_EXTENSION__(),
  applyMiddleware(thunk, updateURLMiddleware)
);

export default store;
