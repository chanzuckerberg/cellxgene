import { createStore, applyMiddleware } from "redux";
import thunk from "redux-thunk";

import cascadeReducers from "./cascade";
import undoable from "./undoable";
import config from "./config";
import userInfo from "./userInfo";
import annoMatrix from "./annoMatrix";
import obsCrossfilter from "./obsCrossfilter";
import categoricalSelection from "./categoricalSelection";
import continuousSelection from "./continuousSelection";
import graphSelection from "./graphSelection";
import colors from "./colors";
import differential from "./differential";
import layoutChoice from "./layoutChoice";
import controls from "./controls";
import annotations from "./annotations";
import genesets from "./genesets";
import genesetsUI from "./genesetsUI";
import autosave from "./autosave";
import ontology from "./ontology";
import centroidLabels from "./centroidLabels";
import pointDialation from "./pointDilation";
import { reembedController } from "./reembed";
import { gcMiddleware as annoMatrixGC } from "../annoMatrix";

import undoableConfig from "./undoableConfig";

const Reducer = undoable(
  cascadeReducers([
    ["config", config],
    ["annoMatrix", annoMatrix],
    ["obsCrossfilter", obsCrossfilter],
    ["ontology", ontology],
    ["annotations", annotations],
    ["genesets", genesets],
    ["genesetsUI", genesetsUI],
    ["layoutChoice", layoutChoice],
    ["categoricalSelection", categoricalSelection],
    ["continuousSelection", continuousSelection],
    ["graphSelection", graphSelection],
    ["colors", colors],
    ["controls", controls],
    ["differential", differential],
    ["centroidLabels", centroidLabels],
    ["pointDilation", pointDialation],
    ["reembedController", reembedController],
    ["autosave", autosave],
    ["userInfo", userInfo],
  ]),
  [
    "annoMatrix",
    "obsCrossfilter",
    "categoricalSelection",
    "continuousSelection",
    "graphSelection",
    "colors",
    "controls",
    "differential",
    "layoutChoice",
    "centroidLabels",
    "annotations",
  ],
  undoableConfig
);

const store = createStore(Reducer, applyMiddleware(thunk, annoMatrixGC));

export default store;
