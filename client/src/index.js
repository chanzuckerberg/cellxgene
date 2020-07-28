import React from "react";
import ReactDOM from "react-dom";
import { Provider } from "react-redux";
import { FocusStyleManager } from "@blueprintjs/core";

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
