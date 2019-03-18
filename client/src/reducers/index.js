// jshint esversion: 6
import { createStore, applyMiddleware } from "redux";
import thunk from "redux-thunk";
import { composeWithDevTools } from "redux-devtools-extension";

import cascadeReducers from "./cascade";
import config from "./config";
import universe from "./universe";
import world from "./world";
import categoricalSelectionState from "./categoricalSelectionState";
import crossfilter from "./crossfilter";
import colors from "./colors";
import differential from "./differential";
import responsive from "./responsive";
import controls from "./controls";
import resetCache from "./resetCache";

const Reducer = cascadeReducers([
  ["config", config],
  ["universe", universe],
  ["world", world],
  ["categoricalSelectionState", categoricalSelectionState],
  ["crossfilter", crossfilter],
  ["colors", colors],
  ["responsive", responsive],
  ["controls", controls],
  ["differential", differential],
  ["resetCache", resetCache]
]);

const store = createStore(Reducer, composeWithDevTools(applyMiddleware(thunk)));

export default store;
