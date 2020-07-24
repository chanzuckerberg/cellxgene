// jshint esversion: 6
const path = require("path");
const webpack = require("webpack");
const HtmlWebpackPlugin = require("html-webpack-plugin");
const FaviconsWebpackPlugin = require("favicons-webpack-plugin");
const ObsoleteWebpackPlugin = require("obsolete-webpack-plugin");

const src = path.resolve("src");
const fonts = path.resolve("src/fonts");
const nodeModules = path.resolve("node_modules");

const babelOptions = require("../babel/babel.dev");

module.exports = {
  mode: "development",
  devtool: "eval",
  entry: ["./src/index"],
  output: {
    path: path.resolve("build"),
    pathinfo: true,
    filename: "static/js/bundle.js",
    publicPath: "/",
  },
  module: {
    rules: [
      {
        test: /\.jsx?$/,
        loader: "babel-loader",
        options: babelOptions,
      },
      {
        test: /\.css$/,
        include: src,
        exclude: [path.resolve(src, "index.css")],
        loader: [
          {
            loader: "style-loader",
          },
          {
            loader: "css-loader",
            options: {
              modules: {
                localIdentName: "[name]__[local]___[hash:base64:5]",
              },
            },
          },
        ],
      },
      {
        test: /index\.css$/,
        include: [path.resolve(src, "index.css")],
        loader: [
          {
            loader: "style-loader",
          },
          {
            loader: "css-loader",
          },
        ],
      },

      { test: /\.json$/, include: [src, nodeModules], loader: "json-loader" },
      {
        test: /\.(jpg|png|gif|eot|svg|ttf|woff|woff2|otf)$/i,
        loader: "file-loader",
        include: [nodeModules, fonts],
        query: { name: "static/assets/[name].[ext]" },
      },
    ],
  },
  plugins: [
    new HtmlWebpackPlugin({
      inject: true,
      template: path.resolve("index.html"),
    }),
    new ObsoleteWebpackPlugin({
      template:
        '<div class="outerContainer">' +
        '<h1 class="purple"> Meta </h1>' +
        '<div class="innerContainer">' +
        '<div class="header"> Unsupported Browser </div>' +
        '<div class="content"> Meta is currently supported on the following browsers' +
        '<div class="supported-browsers purple">' +
        '<div class="browser-type"> <img src="https://assets.beta.meta.org/images/browsers/chrome.png"/>' +
        '<div class="version">Chrome > 60</div> </div>' +
        '<div class="browser-type"> <img src="https://assets.beta.meta.org/images/browsers/safari.png"/>' +
        '<div class="version">Safari >= 10.1</div> </div>' +
        '<div class="browser-type"> <img src="https://assets.beta.meta.org/images/browsers/firefox.png"/>' +
        '<div class="version">Firefox >= 60</div> </div>' +
        '<div class="browser-type"><img src="https://assets.beta.meta.org/images/browsers/edge.png"/>' +
        '<div class="version">Edge >= 15</div> </div>' +
        "</div>",
      promptOnNonTargetBrowser: true,
    }),
    new FaviconsWebpackPlugin({
      logo: "./favicon.png",
      prefix: "static/img/",
      favicons: {
        icons: {
          android: false,
          appleIcon: false,
          appleStartup: false,
          coast: false,
          firefox: false,
          windows: false,
          yandex: false,
        },
      },
    }),
    new webpack.NoEmitOnErrorsPlugin(),
    new webpack.DefinePlugin({
      __REACT_DEVTOOLS_GLOBAL_HOOK__: "({ isDisabled: true })",
    }),
    new webpack.DefinePlugin({
      "process.env.CXG_SERVER_PORT": JSON.stringify(
        process.env.CXG_SERVER_PORT
      ),
    }),
  ],
};
