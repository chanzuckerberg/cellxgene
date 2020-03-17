// jshint esversion: 6
const path = require("path");
const HtmlWebpackPlugin = require("html-webpack-plugin");
const MiniCssExtractPlugin = require("mini-css-extract-plugin");
const HtmlWebpackInlineSourcePlugin = require("html-webpack-inline-source-plugin");
const FaviconsWebpackPlugin = require("favicons-webpack-plugin");
const { CleanWebpackPlugin } = require("clean-webpack-plugin");
const CspHashPlugin = require("./cspHashPlugin");
const SentryCliPlugin = require("@sentry/webpack-plugin");

const src = path.resolve("src");
const fonts = path.resolve("src/fonts");
const nodeModules = path.resolve("node_modules");

const babelOptions = require("../babel/babel.prod");

const publicPath = "/";

const conditionalPlugins = [];

const sentryConfigured = [
  "SENTRY_URL",
  "SENTRY_ORG",
  "SENTRY_PROJECT",
  "SENTRY_AUTH_TOKEN",
].map(varName => varName in process.env).reduce((l, r) => l && r);

if (sentryConfigured) conditionalPlugins.push(
  new SentryCliPlugin({
    include: ".",
    ignoreFile: ".sentrycliignore",
    ignore: ["node_modules", "configuration", "coverage", "__tests__", ".idea"],
    release: `cellxgene@${process.env.CELLXGENE_COMMIT}`
  })
);

module.exports = {
  mode: "production",
  bail: true,
  cache: false,
  entry: ["./src/index.js"],
  output: {
    path: path.resolve("build"),
    publicPath,
  },
  module: {
    rules: [
      {
        test: /\.js$/,
        include: src,
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
        use: [
          MiniCssExtractPlugin.loader,
          {
            loader: "css-loader",
            options: {
              importLoaders: 1,
            },
          },
        ],
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
    ...[
      new HtmlWebpackPlugin({
        inject: "body",
        filename: "index.html",
        template: path.resolve("index_template.html"),
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
          minifyURLs: true,
        },
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
      new HtmlWebpackInlineSourcePlugin(HtmlWebpackPlugin),
      new MiniCssExtractPlugin(),
      new CspHashPlugin({
        filename: "csp-hashes.json",
      }),
    ],
    ...conditionalPlugins
  ],
  performance: {
    maxEntrypointSize: 2000000,
    maxAssetSize: 2000000,
  },
};
