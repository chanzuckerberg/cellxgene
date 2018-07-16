// jshint esversion: 6
/* eslint-disable no-console */
import React from "react";
import ReactDOM from "react-dom";
import { AppContainer } from "react-hot-loader";
import { Provider } from "react-redux";
import Redbox from "redbox-react";

/* our code */
import App from "./components/app";
import store from "./reducers";

ReactDOM.render(
  <AppContainer errorReporter={Redbox}>
    <Provider store={store}>
      <App />
    </Provider>
  </AppContainer>,
  document.getElementById("root")
);

// Hot Module Replacement API
if (module.hot) {
  module.hot.accept("./components/app", () => {
    const NextApp = require("./components/app").default;
    ReactDOM.render(
      <AppContainer>
        <Provider store={store}>
          <NextApp />
        </Provider>
      </AppContainer>,
      document.getElementById("root")
    );
  });
}
