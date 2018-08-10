// jshint esversion: 6
const path = require("path");
const autoprefixer = require("autoprefixer");
const webpack = require("webpack");
const HtmlWebpackPlugin = require("html-webpack-plugin");
const MiniCssExtractPlugin = require("mini-css-extract-plugin");
const SWPrecacheWebpackPlugin = require("sw-precache-webpack-plugin");
const CopyWebpackPlugin = require("copy-webpack-plugin");
const HtmlWebpackInlineSourcePlugin = require("html-webpack-inline-source-plugin");
const MinifyPlugin = require("babel-minify-webpack-plugin");

const src = path.resolve("src");
const nodeModules = path.resolve("node_modules");

const publicPath = "/";

module.exports = {
  mode: "production",
  bail: true,
  // TODO: causes a js error, need to update to whatever webpack4 wants
  // devtool: "cheap-source-map",
  entry: [require.resolve("../polyfills/polyfills"), path.join(src, "index")],
  output: {
    path: path.resolve("build"),
    filename: "static/js/[name].[chunkhash:8].js",
    chunkFilename: "static/js/[name].[chunkhash:8].chunk.js",
    publicPath
  },
  resolve: { extensions: [".js", ".json"] },
  module: {
    rules: [
      {
        test: /\.js$/,
        include: src,
        loader: "babel-loader",
        query: require("../babel/babel.prod")
      },
      {
        test: /\.css$/,
        include: [src, nodeModules],
        loader: [
          {
            loader: MiniCssExtractPlugin.loader
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
              plugins: function() {
                return [];
              }
            }
          }
        ]
      },
      {
        test: /\.json$/,
        include: [src, nodeModules],
        loader: "json-loader",
        exclude: /manifest.json$/
      },
      {
        test: /\.(jpg|png|gif|eot|svg|ttf|woff|woff2)(\?.*)?$/,
        include: [src, nodeModules],
        loader: "file-loader",
        query: { name: "static/media/[name].[hash:8].[ext]" }
      },
      {
        test: /\.(mp4|webm)(\?.*)?$/,
        include: [src, nodeModules],
        loader: "url-loader",
        query: { limit: 10000, name: "static/media/[name].[hash:8].[ext]" }
      }
    ]
  },
  plugins: [
    new HtmlWebpackPlugin({
      inject: "body",
      filename: "index.html",
      template: path.resolve("index_template.html"),
      favicon: path.resolve("favicon.png"),
      inlineSource: ".(js|css)$",
      minify: {
        removeComments: true,
        collapseWhitespace: true,
        removeRedundantAttributes: true,
        useShortDoctype: true,
        removeEmptyAttributes: true,
        removeStyleLinkTypeAttributes: true,
        keepClosingSlash: true,
        minifyJS: true,
        minifyCSS: true,
        minifyURLs: true
      }
    }),
    new HtmlWebpackInlineSourcePlugin(),
    new webpack.LoaderOptionsPlugin({
      options: {
        eslint: {
          configFile: path.resolve("./configuration/eslint/eslint.js"),
          useEslintrc: false
        },
        postcss() {
          return [autoprefixer];
        }
      }
    }),
    new webpack.DefinePlugin({ "process.env.NODE_ENV": '"production"' }),
    new webpack.optimize.OccurrenceOrderPlugin(),
    new MinifyPlugin(),
    new MiniCssExtractPlugin({
      filename: "static/css/[name].[contenthash:8].css"
    }),
    new SWPrecacheWebpackPlugin({
      cacheId: "cellxgene",
      filename: "service-worker.js"
    })
  ]
};
