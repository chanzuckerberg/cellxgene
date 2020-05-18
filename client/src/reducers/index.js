import { createStore, applyMiddleware } from "redux";
import thunk from "redux-thunk";

import cascadeReducers from "./cascade";
import undoable from "./undoable";
import config from "./config";
import universe from "./universe";
import world from "./world";
import categoricalSelection from "./categoricalSelection";
import continuousSelection from "./continuousSelection";
import graphSelection from "./graphSelection";
import crossfilter from "./crossfilter";
import colors from "./colors";
import differential from "./differential";
import layoutChoice from "./layoutChoice";
import controls from "./controls";
import resetCache from "./resetCache";
import annotations from "./annotations";
import autosave from "./autosave";
import ontology from "./ontology";
import centroidLabels from "./centroidLabels";
import pointDialation from "./pointDilation";
import { reembedController, reembedding } from "./reembed";

import undoableConfig from "./undoableConfig";

const Reducer = undoable(
  cascadeReducers([
    ["config", config],
    ["universe", universe],
    ["world", world],
    ["ontology", ontology],
    ["annotations", annotations],
    ["reembedding", reembedding],
    ["layoutChoice", layoutChoice],
    ["categoricalSelection", categoricalSelection],
    ["continuousSelection", continuousSelection],
    ["graphSelection", graphSelection],
    ["colors", colors],
    ["controls", controls],
    ["crossfilter", crossfilter],
    ["differential", differential],
    ["centroidLabels", centroidLabels],
    ["pointDilation", pointDialation],
    ["reembedController", reembedController],
    ["autosave", autosave],
    ["resetCache", resetCache],
  ]),
  [
    "universe",
    "world",
    "categoricalSelection",
    "continuousSelection",
    "graphSelection",
    "crossfilter",
    "colors",
    "controls",
    "differential",
    "layoutChoice",
    "centroidLabels",
    "annotations",
    "reembedding",
  ],
  undoableConfig
);

const store = createStore(Reducer, applyMiddleware(thunk));

export default store;
