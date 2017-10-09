import { combineReducers, createStore, applyMiddleware } from 'redux';
import updateURLMiddleware from "../middleware/updateURLMiddleware";
import cells from './cells';
import onURLchange from "./onURLchange";

const Reducer = combineReducers({
  cells,
  onURLchange,
})

let store = createStore(
  Reducer,
  window.__REDUX_DEVTOOLS_EXTENSION__ && window.__REDUX_DEVTOOLS_EXTENSION__(),
  applyMiddleware(updateURLMiddleware)
);

export default store;
