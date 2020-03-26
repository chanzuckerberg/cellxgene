// jshint esversion: 6
const path = require("path");
const webpack = require("webpack");
const HtmlWebpackPlugin = require("html-webpack-plugin");

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
    publicPath: "/"
  },
  module: {
    rules: [
      {
        test: /\.js$/,
        include: src,
        loader: "babel-loader",
        options: babelOptions
      },
      {
        test: /\.css$/,
        include: src,
        exclude: [path.resolve(src, "index.css")],
        loader: [
          {
            loader: "style-loader"
          },
          {
            loader: "css-loader",
            options: {
              modules: {
                localIdentName: "[name]__[local]___[hash:base64:5]"
              }
            }
          }
        ]
      },
      {
        test: /index\.css$/,
        include: [path.resolve(src, "index.css")],
        loader: [
          {
            loader: "style-loader"
          },
          {
            loader: "css-loader"
          }
        ]
      },

      { test: /\.json$/, include: [src, nodeModules], loader: "json-loader" },
      {
        test: /\.(jpg|png|gif|eot|svg|ttf|woff|woff2)(\?.*)?$/,
        include: nodeModules,
        loader: "file-loader",
        query: { name: "static/media/[name].[ext]" }
      },
      {
        test: /\.(woff|woff2|eot|ttf|otf)$/,
        include: fonts,
        loader: "file-loader",
        query: { name: "static/fonts/[name].[ext]" }
      },
      {
        test: /\.(mp4|webm)(\?.*)?$/,
        include: [src, nodeModules],
        loader: "url-loader",
        query: { limit: 10000, name: "static/media/[name].[ext]" }
      }
    ]
  },
  plugins: [
    new HtmlWebpackPlugin({
      inject: true,
      template: path.resolve("index.html"),
      favicon: path.resolve("favicon.png")
    }),
    new webpack.NoEmitOnErrorsPlugin(),
    new webpack.DefinePlugin({
      __REACT_DEVTOOLS_GLOBAL_HOOK__: "({ isDisabled: true })"
    }),
    new webpack.DefinePlugin({
      "process.env.CXG_SERVER_PORT": JSON.stringify(process.env.CXG_SERVER_PORT)
    })
  ]
};
