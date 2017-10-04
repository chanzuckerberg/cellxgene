import { combineReducers, createStore } from 'redux'
import cells from './cells'

const Reducer = combineReducers({
  cells,
})

let store = createStore(Reducer);

export default store;
