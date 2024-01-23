const path = require("path");
const webpack = require("webpack");
const HtmlWebpackPlugin = require("html-webpack-plugin");
const { CleanWebpackPlugin } = require("clean-webpack-plugin");
const TerserJSPlugin = require("terser-webpack-plugin");
const FaviconsWebpackPlugin = require("favicons-webpack-plugin");
const MiniCssExtractPlugin = require("mini-css-extract-plugin");

const { merge } = require("webpack-merge");

const babelOptions = require("../babel/babel.prod");

const CspHashPlugin = require("./cspHashPlugin");
const sharedConfig = require("./webpack.config.shared");

const fonts = path.resolve("src/fonts");
const nodeModules = path.resolve("node_modules");

const prodConfig = {
  mode: "production",
  bail: true,
  cache: false,
  output: {
    filename: "static/[name]-[contenthash].js",
  },
  optimization: {
    minimize: true,
    minimizer: [new TerserJSPlugin({})],
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
        test: /\.(jpg|png|gif|eot|svg|ttf|woff|woff2|otf)$/i,
        loader: "file-loader",
        include: [nodeModules, fonts],
        options: {
          name: "static/assets/[name]-[contenthash].[ext]",
          // (thuang): This is needed to make sure @font url path is '../static/assets/'
          publicPath: "static/",
        },
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
    new webpack.DefinePlugin({
      // webpack 5 no longer polyfills NodeJS modules, so fake the one we need
      "process.env": JSON.stringify({
        NODE_ENV: "production",
      }),
    }),
  ],
  performance: {
    maxEntrypointSize: 2000000,
    maxAssetSize: 2000000,
  },
};

module.exports = merge(sharedConfig, prodConfig);
