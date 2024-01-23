const path = require("path");
const webpack = require("webpack");
const HtmlWebpackPlugin = require("html-webpack-plugin");
const FaviconsWebpackPlugin = require("favicons-webpack-plugin");
const MiniCssExtractPlugin = require("mini-css-extract-plugin");

const { merge } = require("webpack-merge");

const sharedConfig = require("./webpack.config.shared");
const babelOptions = require("../babel/babel.dev");

const fonts = path.resolve("src/fonts");
const nodeModules = path.resolve("node_modules");

const devConfig = {
  mode: "development",
  devtool: "eval",
  output: {
    pathinfo: true,
    filename: "static/js/bundle.js",
  },
  module: {
    rules: [
      {
        test: /\.jsx?$/,
        loader: "babel-loader",
        options: babelOptions,
      },
      {
        test: /\.(jpg|png|gif|eot|svg|ttf|woff|woff2|otf)$/i,
        loader: "file-loader",
        include: [nodeModules, fonts],
        options: {
          name: "static/assets/[name].[ext]",
          // (thuang): This is needed to make sure @font url path is '/static/assets/'
          publicPath: "/",
        },
      },
    ],
  },
  plugins: [
    new HtmlWebpackPlugin({
      inject: true,
      template: path.resolve("index.html"),
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
    new MiniCssExtractPlugin({
      filename: "static/[name].css",
    }),
    new webpack.NoEmitOnErrorsPlugin(),
    new webpack.DefinePlugin({
      __REACT_DEVTOOLS_GLOBAL_HOOK__: "({ isDisabled: true })",
    }),
    new webpack.DefinePlugin({
      // webpack 5 no longer polyfills NodeJS modules, so fake the one we need
      "process.env": JSON.stringify({
        NODE_ENV: process.env.NODE_ENV || "development",
        CXG_SERVER_PORT: process.env.CXG_SERVER_PORT || "5005",
      }),
    }),
  ],
  infrastructureLogging: {
    level: "warn",
  },
};

module.exports = merge(sharedConfig, devConfig);
