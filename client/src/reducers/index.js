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
import responsive from "./responsive";
import controls from "./controls";
import resetCache from "./resetCache";
import annotations from "./annotations";
import autosave from "./autosave";
import centroidLabels from "./centroidLabels";
import pointDialation from "./pointDilation";

import undoableConfig from "./undoableConfig";

const Reducer = undoable(
    cascadeReducers([
        ["universe", universe],
        ["config", config],
        ["world", world],
        ["annotations", annotations],
        ["layoutChoice", layoutChoice],
        ["categoricalSelection", categoricalSelection],
        ["continuousSelection", continuousSelection],
        ["graphSelection", graphSelection],
        ["crossfilter", crossfilter],
        ["colors", colors],
        ["controls", controls],
        ["differential", differential],
        ["responsive", responsive],
    ["centroidLabels", centroidLabels],
    ["pointDilation", pointDialation],
        ["autosave", autosave],
        ["resetCache", resetCache]
    ]),
    [
        "universe",
        "categoricalSelection",
        "world",
        "continuousSelection",
        "graphSelection",
        "crossfilter",
        "controls",
        "differential",
        "layoutChoice",
        "colors",
        "annotations"
);

const store = createStore(Reducer, applyMiddleware(thunk));

export default store;
