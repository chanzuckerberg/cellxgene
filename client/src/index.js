// jshint esversion: 6
/* eslint-disable no-console */
import React from "react";
import ReactDOM from "react-dom";
import { Provider } from "react-redux";

/* our code */
import App from "./components/app";
import store from "./reducers";

ReactDOM.render(
  <Provider store={store}>
    <App />
  </Provider>,
  document.getElementById("root")
);
