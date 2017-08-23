/* eslint-disable */
var path = require('path');
var historyApiFallback = require('connect-history-api-fallback');
var chalk = require('chalk');
var express = require('express');
var favicon = require('serve-favicon');
var webpack = require('webpack');
var config = require('../configuration/webpack/webpack.config.dev');
var utils = require('./utils');

process.env.NODE_ENV = 'development';

var PORT = process.env.PORT || 3000;

// Set up compiler
var compiler = webpack(config);

compiler.plugin('invalid', () => {
  utils.clearConsole();
  console.log('Compiling...');
});

compiler.plugin('done', stats => {
  utils.formatStats(stats, PORT);
});

// Launch server
var app = express();

app.use(historyApiFallback({ verbose: false }));

app.use(
  require('webpack-dev-middleware')(compiler, {
    noInfo: true,
    publicPath: config.output.publicPath
  })
);

app.use(require('webpack-hot-middleware')(compiler));

app.use(favicon('./favicon.png'));

app.get('*', (req, res) => {
  res.sendFile(path.resolve('index.html'));
});

app.listen(PORT, err => {
  if (err) {
    console.log(err);
    return;
  }

  utils.clearConsole();
  console.log(chalk.cyan('Starting the development server...'));
  console.log();
});
