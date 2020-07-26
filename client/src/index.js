// jshint esversion: 6
import "core-js/stable";
import "regenerator-runtime/runtime";
import React from "react";
import ReactDOM from "react-dom";
import { Provider } from "react-redux";
import { FocusStyleManager } from "@blueprintjs/core";
// https://github.com/anonyco/FastestSmallestTextEncoderDecoder
// Only supports UTF-8
import "fastestsmallesttextencoderdecoder";

import "./index.css";

/* our code */
import App from "./components/app";
import store from "./reducers";

FocusStyleManager.onlyShowFocusOnTabs();

ReactDOM.render(
  <Provider store={store}>
    <App />
  </Provider>,
  document.getElementById("root")
);
