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
const ScriptExtHtmlWebpackPlugin = require("script-ext-html-webpack-plugin");

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
        "<div>" +
        '<div style="display: flex;flex-direction: column;width: 100vw;height: 100vh;text-align: center;justify-content: center;">' +
        "<h1> cellxgene </h1>" +
        "<div> Unsupported Browser </div>" +
        '<div style="margin-top: 16px;"> cellxgene is currently supported on the following browsers' +
        '<div style="display: flex;justify-content: space-around;margin: auto 40vw;margin-top: 16px">' +
        '<div> <img src="https://assets.beta.meta.org/images/browsers/chrome.png" style="width: 80px;height: 80px;"><div>Chrome &gt; 60</div></div>' +
        '<div> <img src="https://assets.beta.meta.org/images/browsers/safari.png" style="width: 80px;height: 80px;"><div>Safari &gt;= 10.1</div> </div>' +
        '<div> <img src="https://assets.beta.meta.org/images/browsers/firefox.png" style="width: 80px;height: 80px;"><div>Firefox &gt;= 60</div> </div>' +
        '<div><img src="https://assets.beta.meta.org/images/browsers/edge.png" style="width: 80px;height: 80px;"><div>Edge &gt;= 15</div> </div>' +
        "</div></div><div></div></div>",
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
    new ScriptExtHtmlWebpackPlugin({
      async: "obsolete",
    }),
  ],
  performance: {
    maxEntrypointSize: 2000000,
    maxAssetSize: 2000000,
  },
};
