// jshint esversion: 6
const path = require("path");
const webpack = require("webpack");
const HtmlWebpackPlugin = require("html-webpack-plugin");
const FaviconsWebpackPlugin = require("favicons-webpack-plugin");
const ObsoleteWebpackPlugin = require("obsolete-webpack-plugin");
const ScriptExtHtmlWebpackPlugin = require("script-ext-html-webpack-plugin");

const src = path.resolve("src");
const fonts = path.resolve("src/fonts");
const nodeModules = path.resolve("node_modules");

const babelOptions = require("../babel/babel.dev");

module.exports = {
  mode: "development",
  devtool: "eval",
  entry: ["whatwg-fetch", "abort-controller/polyfill", "./src/index"],
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
      name: "obsolete",
      template:
        "<script>" +
        'var root = document.getElementById("root");' +
        "root.remove();" +
        'var portals = document.getElementsByClassName("bp3-portal");' +
        "for(var i = 0; i < portals.length; i += 1) {" +
        "portals[i].remove();" +
        "}" +
        "</script>" +
        '<div style="display: flex;flex-direction: column;width: 100vw;height: 100vh;text-align: center;justify-content: center;align-items: center;background: #8080801a;font-family: &quot;Roboto Condensed, sans serif&quot;;">' +
        '<img src="https://raw.githubusercontent.com/chanzuckerberg/cellxgene/main/docs/cellxgene-logo.png" style="width: 320px;"/>' +
        '<div style="margin-top: 16px;background: white;width: 40vw;border-radius: 4px;padding: 24px 64px;-webkit-box-shadow: 0px 0px 3px 2px rgba(0,0,0,0.38);-moz-box-shadow: 0px 0px 3px 2px rgba(0,0,0,0.38);box-shadow: 0px 0px 3px 2px rgba(0,0,0,0.38);max-width: 550px;">' +
        '<div style="margin-bottom: 0;font-weight: bolder;font-size: 1.2em;"> Unsupported Browser </div>' +
        '<div style="margin-top: 0;">cellxgene is currently supported on the following browsers</div>' +
        '<div style="display: flex;justify-content: space-around;margin-top: 16px">' +
        "<div>" +
        '<img src="https://assets.beta.meta.org/images/browsers/chrome.png" style="width: 80px;height: 80px;"/>' +
        "<div>Chrome &gt; 60</div>" +
        "</div>" +
        "<div>" +
        '<img src="https://assets.beta.meta.org/images/browsers/safari.png" style="width: 80px;height: 80px;"/>' +
        "<div>Safari ≥ 10.1</div>" +
        "</div>" +
        "<div>" +
        '<img src="https://assets.beta.meta.org/images/browsers/firefox.png" style="width: 80px;height: 80px;"/>' +
        "<div>Firefox ≥ 60</div>" +
        "</div>" +
        "<div>" +
        '<img src="https://assets.beta.meta.org/images/browsers/edge.png" style="width: 80px;height: 80px;"/>' +
        "<div>Edge ≥ 15</div>" +
        "</div>" +
        "</div>" +
        "</div>" +
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
    new ScriptExtHtmlWebpackPlugin({
      async: "obsolete",
    }),
  ],
};
