const path = require('path');
const autoprefixer = require('autoprefixer');
const webpack = require('webpack');
const HtmlWebpackPlugin = require('html-webpack-plugin');
const ExtractTextPlugin = require('extract-text-webpack-plugin');
const SWPrecacheWebpackPlugin = require('sw-precache-webpack-plugin');
const CopyWebpackPlugin = require('copy-webpack-plugin');
const HtmlWebpackInlineSourcePlugin = require(
  'html-webpack-inline-source-plugin'
);

const src = path.resolve('src');
const nodeModules = path.resolve('node_modules');

const publicPath = '/';

module.exports = {
  bail: true,
  devtool: 'source-map',
  entry: [ require.resolve('../polyfills/polyfills'), path.join(src, 'index') ],
  output: {
    path: path.resolve('build'),
    filename: 'static/js/[name].[chunkhash:8].js',
    chunkFilename: 'static/js/[name].[chunkhash:8].chunk.js',
    publicPath
  },
  resolve: { extensions: [ '.js', '.json' ] },
  module: {
    loaders: [
      {
        test: /\.js$/,
        include: src,
        loader: 'babel-loader',
        query: require('../babel/babel.prod')
      },
      {
        test: /\.css$/,
        include: [ src, nodeModules ],
        loader: ExtractTextPlugin.extract({
          fallbackLoader: 'style-loader',
          loader: 'css-loader?modules&importLoaders=1&localIdentName=[name]__[local]___[hash:base64:5]-autoprefixer!postcss-loader'
        })
      },
      {
        test: /\.json$/,
        include: [ src, nodeModules ],
        loader: 'json-loader',
        exclude: /manifest.json$/
      },
      {
        test: /\.(jpg|png|gif|eot|svg|ttf|woff|woff2)(\?.*)?$/,
        include: [ src, nodeModules ],
        loader: 'file-loader',
        query: { name: 'static/media/[name].[hash:8].[ext]' }
      },
      {
        test: /\.(mp4|webm)(\?.*)?$/,
        include: [ src, nodeModules ],
        loader: 'url-loader',
        query: { limit: 10000, name: 'static/media/[name].[hash:8].[ext]' }
      }
    ]
  },
  plugins: [
    new HtmlWebpackPlugin({
      inject: 'body',
      template: path.resolve('index.html'),
      favicon: path.resolve('favicon.png'),
      inlineSource: '.(js|css)$',
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
          configFile: path.resolve('./configuration/eslint/eslint.js'),
          useEslintrc: false
        },
        postcss() {
          return [ autoprefixer ];
        }
      }
    }),
    new webpack.DefinePlugin({ 'process.env.NODE_ENV': '"production"' }),
    new webpack.optimize.OccurrenceOrderPlugin(),
    new webpack.optimize.UglifyJsPlugin({
      compress: { screw_ie8: true, warnings: false },
      mangle: { screw_ie8: true },
      output: { comments: false, screw_ie8: true }
    }),
    new ExtractTextPlugin('static/css/[name].[contenthash:8].css'),
    new CopyWebpackPlugin([
      { from: 'public' },
      { from: 'manifest.webmanifest' }
    ]),
    new SWPrecacheWebpackPlugin({
      cacheId: 'formidable-react-starter',
      filename: 'service-worker.js'
    })
  ]
};
