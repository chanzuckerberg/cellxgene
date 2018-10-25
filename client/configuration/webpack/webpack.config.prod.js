// jshint esversion: 6
const path = require("path");
const HtmlWebpackPlugin = require("html-webpack-plugin");
const MiniCssExtractPlugin = require("mini-css-extract-plugin");
const SWPrecacheWebpackPlugin = require("sw-precache-webpack-plugin");
const HtmlWebpackInlineSourcePlugin = require("html-webpack-inline-source-plugin");

const src = path.resolve("src");
const nodeModules = path.resolve("node_modules");

const babelOptions = require("../babel/babel.prod");

const publicPath = "/";

module.exports = {
  mode: "production",
  bail: true,
  cache: false,
  devtool: "cheap-source-map",
  entry: ["./src/index.js"],
  output: {
    path: path.resolve("build"),
    filename: "static/js/[name].[chunkhash:8].js",
    publicPath
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
        use: [
          MiniCssExtractPlugin.loader,
          {
            loader: "css-loader",
            options: {
              modules: true,
              importLoaders: 1,
              localIdentName: "[name]__[local]___[hash:base64:5]"
            }
          }
        ]
      },
      {
        test: /index\.css$/,
        include: [path.resolve(src, "index.css")],
        use: [
          MiniCssExtractPlugin.loader,
          {
            loader: "css-loader",
            options: {
              importLoaders: 1
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
        options: {
          name: "[name].[hash:8].[ext]",
          outputPath: "static/media/"
        }
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
    new MiniCssExtractPlugin({
      filename: "static/css/[name].[contenthash:8].css"
    }),
    new SWPrecacheWebpackPlugin({
      cacheId: "cellxgene",
      filename: "service-worker.js"
    })
  ],
  performance: {
    maxEntrypointSize: 2000000,
    maxAssetSize: 2000000
  }
};
