// jshint esversion: 6
import React from "react";
import ReactDOM from "react-dom";
import { Provider } from "react-redux";
import { FocusStyleManager } from "@blueprintjs/core";

import "./index.css";

/* our code */
import App from "./components/app";
import store from "./reducers";

// import * as serviceWorker from "./serviceWorker";
import { Auth0Provider } from "./react-auth0-spa";
// import config from "./auth_config.json";
import history from "./util/history";

// A function that routes the user to the right place
// after login
const onRedirectCallback = (appState) => {
  history.push(
    appState && appState.targetUrl
      ? appState.targetUrl
      : window.location.pathname
  );
};
FocusStyleManager.onlyShowFocusOnTabs();

ReactDOM.render(
  <Auth0Provider
    domain="czi-single-cell.auth0.com"
    client_id="e9VV5VPAa5AVaHRGaFQIJOhc2N2jXNZ6"
    redirect_uri={window.location.origin}
    onRedirectCallback={onRedirectCallback}
  >
    <Provider store={store}>
      <App />
    </Provider>
  </Auth0Provider>,

  document.getElementById("root")
);

// serviceWorker.unregister();
