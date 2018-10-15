// jshint esversion: 6
const path = require("path");
const autoprefixer = require("autoprefixer");
const webpack = require("webpack");
const HtmlWebpackPlugin = require("html-webpack-plugin");

const src = path.resolve("src");
const nodeModules = path.resolve("node_modules");

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
        options: require("../babel/babel.dev")
      },
      {
        test: /\.css$/,
        include: [src, nodeModules],
        loader: [
          {
            loader: "style-loader"
          },
          {
            loader: "css-loader",
            options: {
              modules: true,
              importLoaders: 1,
              localIdentName: "[name]__[local]___[hash:base64:5]"
            }
          },
          {
            loader: "postcss-loader",
            options: {
              plugins: () => []
            }
          }
        ]
      },
      { test: /\.json$/, include: [src, nodeModules], loader: "json-loader" },
      {
        test: /\.(jpg|png|gif|eot|svg|ttf|woff|woff2)(\?.*)?$/,
        include: [src, nodeModules],
        loader: "file-loader",
        query: { name: "static/media/[name].[ext]" }
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
    new webpack.NoEmitOnErrorsPlugin()
  ]
};
