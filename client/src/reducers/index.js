import { createStore, applyMiddleware } from "redux";
import thunk from "redux-thunk";

import cascadeReducers from "./cascade";
import undoable from "./undoable";
import config from "./config";
import universe from "./universe";
import world from "./world";
import categoricalSelection from "./categoricalSelection";
import crossfilter from "./crossfilter";
import colors from "./colors";
import differential from "./differential";
import responsive from "./responsive";
import controls from "./controls";
import resetCache from "./resetCache";

const ignoredActions = new Set([
  // these actions will not affect history, ie, we will
  // not snapshot history upon these actions.  These take
  // precedent over `clearHistoryUponActions`
  "url changed",
  "interface reset started",
  "initial data load start",
  "configuration load complete",
  "increment graph render counter",
  "window resize",

  "lasso started",

  "request differential expression success",

  "expression load start",
  "expression load success",
  "expression load error",

  "continuous metadata histogram brush",
  "continuous metadata histogram end",

  "request user defined gene started",
  "request user defined gene success",
  "request user defined gene error",
  "bulk user defined gene complete",
  "single user defined gene complete"
]);

const clearOnActions = new Set([
  // history will be cleared when these actions occur
  "initial data load complete (universe exists)",
  "reset World to eq Universe",
  "initial data load error"
]);

/* configuration for the undoable meta reducer */
const undoableConfig = {
  historyLimit: 50, // maximum history size
  skipActionFilter: (state, action) => ignoredActions.has(action.type),
  clearOnActionFilter: (state, action) => clearOnActions.has(action.type)
};

const Reducer = undoable(
  cascadeReducers([
    ["config", config],
    ["universe", universe],
    ["world", world],
    ["categoricalSelection", categoricalSelection],
    ["crossfilter", crossfilter],
    ["colors", colors],
    ["controls", controls],
    ["differential", differential],
    ["responsive", responsive],
    ["resetCache", resetCache]
  ]),
  [
    "world",
    "categoricalSelection",
    "crossfilter",
    "colors",
    "controls",
    "differential"
  ],
  undoableConfig
);

const store = createStore(Reducer, applyMiddleware(thunk));

export default store;
