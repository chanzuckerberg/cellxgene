/* eslint-disable no-console */
import React from 'react';
import ReactDOM from 'react-dom';
import { AppContainer } from 'react-hot-loader';
import Redbox from 'redbox-react';

import Routes from './routes';
import './index.css';

if (process.env.NODE_ENV !== 'development') {
  if ('serviceWorker' in navigator) {
    navigator.serviceWorker.register('/service-worker.js');
  }

  ReactDOM.render(<Routes />, document.getElementById('root'));
} else {
  ReactDOM.render(
    <AppContainer errorReporter={Redbox}>
      <Routes />
    </AppContainer>,
    document.getElementById('root')
  );

  if (module.hot) {
    module.hot.accept('./routes', () => {
      // TODO: Remove console override when
      // https://github.com/reactjs/react-router/issues/2704 is fixed
      const orgError = console.error;
      console.error = message => {
        if (
          message &&
            message.indexOf('You cannot change <Router routes>;') === -1
        ) {
          orgError.apply(console, [ message ]);
        }
      };

      const RoutesUpdate = require('./routes').default;

      ReactDOM.render(
        <AppContainer errorReporter={Redbox}>
          <RoutesUpdate />
        </AppContainer>,
        document.getElementById('root')
      );
    });
  }
}
