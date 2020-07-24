// jshint esversion: 6
const path = require("path");
const HtmlWebpackPlugin = require("html-webpack-plugin");
const MiniCssExtractPlugin = require("mini-css-extract-plugin");
const FaviconsWebpackPlugin = require("favicons-webpack-plugin");
const { CleanWebpackPlugin } = require("clean-webpack-plugin");
const TerserJSPlugin = require("terser-webpack-plugin");
const CleanCss = require("clean-css");
const OptimizeCSSAssetsPlugin = require("optimize-css-assets-webpack-plugin");
const ObsoleteWebpackPlugin = require("obsolete-webpack-plugin");

const CspHashPlugin = require("./cspHashPlugin");

const src = path.resolve("src");
const fonts = path.resolve("src/fonts");
const nodeModules = path.resolve("node_modules");

const babelOptions = require("../babel/babel.prod");

const publicPath = "/";

module.exports = {
  mode: "production",
  bail: true,
  cache: false,
  entry: ["./src/index.js"],
  output: {
    filename: "static/[name]-[contenthash].js",
    path: path.resolve("build"),
    publicPath,
  },
  optimization: {
    minimize: true,
    minimizer: [
      new TerserJSPlugin({}),
      new OptimizeCSSAssetsPlugin({
        cssProcessor: CleanCss,
      }),
    ],
  },
  devtool: "source-map",
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
        use: [
          MiniCssExtractPlugin.loader,
          {
            loader: "css-loader",
            options: {
              modules: {
                localIdentName: "[name]__[local]___[hash:base64:5]",
              },
              importLoaders: 1,
            },
          },
        ],
      },
      {
        test: /index\.css$/,
        include: [path.resolve(src, "index.css")],
        use: [MiniCssExtractPlugin.loader, "css-loader"],
      },
      {
        test: /\.json$/,
        include: [src, nodeModules],
        loader: "json-loader",
        exclude: /manifest.json$/,
      },
      {
        test: /\.(jpg|png|gif|eot|svg|ttf|woff|woff2|otf)$/i,
        loader: "file-loader",
        include: [nodeModules, fonts],
        query: { name: "static/assets/[name]-[contenthash].[ext]" },
      },
    ],
  },
  plugins: [
    new HtmlWebpackPlugin({
      filename: "index.html",
      template: path.resolve("index_template.html"),
      decodeEntities: false,
      minify: false,
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
    new CleanWebpackPlugin({
      verbose: true,
      protectWebpackAssets: false,
      cleanAfterEveryBuildPatterns: ["main.js", "main.css"],
    }),
    new FaviconsWebpackPlugin({
      logo: "./favicon.png",
      prefix: "static/assets/",
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
      filename: "static/[name]-[contenthash].css",
    }),
    new CspHashPlugin({
      filename: "csp-hashes.json",
    }),
  ],
  performance: {
    maxEntrypointSize: 2000000,
    maxAssetSize: 2000000,
  },
};
