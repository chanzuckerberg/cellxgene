/* eslint-disable @blueprintjs/classes-constants -- we don't import blueprint here  */
// eslint-disable-next-line @typescript-eslint/no-var-requires --- FIXME: disabled temporarily on migrate to TS.
const path = require("path");
// eslint-disable-next-line @typescript-eslint/no-var-requires --- FIXME: disabled temporarily on migrate to TS.
const fs = require("fs");
// eslint-disable-next-line @typescript-eslint/no-var-requires --- FIXME: disabled temporarily on migrate to TS.
const MiniCssExtractPlugin = require("mini-css-extract-plugin");
// eslint-disable-next-line @typescript-eslint/no-var-requires --- FIXME: disabled temporarily on migrate to TS.
const ObsoleteWebpackPlugin = require("obsolete-webpack-plugin");
// eslint-disable-next-line @typescript-eslint/no-var-requires --- FIXME: disabled temporarily on migrate to TS.
const ScriptExtHtmlWebpackPlugin = require("script-ext-html-webpack-plugin");

const src = path.resolve("src");
const nodeModules = path.resolve("node_modules");

const publicPath = "";

const rawObsoleteHTMLTemplate = fs.readFileSync(
  `${__dirname}/obsoleteHTMLTemplate.html`,
  "utf8"
);

const obsoleteHTMLTemplate = rawObsoleteHTMLTemplate.replace(/'/g, '"');

module.exports = {
  entry: [
    "core-js",
    "regenerator-runtime/runtime",
    "whatwg-fetch",
    "abort-controller/polyfill",
    "./src/index",
  ],
  output: {
    path: path.resolve("build"),
    publicPath,
  },
  resolve: {
    extensions: [".ts", ".tsx", "..."],
  },
  module: {
    rules: [
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
                localIdentName: "[name]__[local]___[contenthash:base64:5]",
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
    ],
  },
  plugins: [
    new ObsoleteWebpackPlugin({
      name: "obsolete",
      template: obsoleteHTMLTemplate,
      promptOnNonTargetBrowser: false,
    }),
    new ScriptExtHtmlWebpackPlugin({
      async: "obsolete",
    }),
  ],
};
/* eslint-enable @blueprintjs/classes-constants -- we don't import blueprint here  */
