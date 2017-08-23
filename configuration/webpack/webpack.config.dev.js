const path = require('path');
const autoprefixer = require('autoprefixer');
const webpack = require('webpack');
const HtmlWebpackPlugin = require('html-webpack-plugin');

const src = path.resolve('src');
const nodeModules = path.resolve('node_modules');

module.exports = {
  devtool: 'eval',
  entry: [
    'webpack-hot-middleware/client?quiet=true&noInfo=true',
    require.resolve('react-hot-loader/patch'),
    require.resolve('../polyfills/polyfills'),
    path.join(src, 'index')
  ],
  output: {
    path: path.resolve('build'),
    pathinfo: true,
    filename: 'static/js/bundle.js',
    publicPath: '/'
  },
  resolve: { extensions: [ '.js', '.json' ] },
  module: {
    loaders: [
      {
        test: /\.js$/,
        include: src,
        loader: 'babel-loader',
        query: require('../babel/babel.dev')
      },
      {
        test: /\.css$/,
        include: [ src, nodeModules ],
        loader: 'style-loader!css-loader?modules&importLoaders=1&localIdentName=[name]__[local]___[hash:base64:5]!postcss-loader'
      },
      { test: /\.json$/, include: [ src, nodeModules ], loader: 'json-loader' },
      {
        test: /\.(jpg|png|gif|eot|svg|ttf|woff|woff2)(\?.*)?$/,
        include: [ src, nodeModules ],
        loader: 'file-loader',
        query: { name: 'static/media/[name].[ext]' }
      },
      {
        test: /\.(mp4|webm)(\?.*)?$/,
        include: [ src, nodeModules ],
        loader: 'url-loader',
        query: { limit: 10000, name: 'static/media/[name].[ext]' }
      }
    ]
  },
  plugins: [
    new HtmlWebpackPlugin({
      inject: true,
      template: path.resolve('index.html'),
      favicon: path.resolve('favicon.png')
    }),
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
    new webpack.DefinePlugin({ 'process.env.NODE_ENV': '"development"' }),
    // Note: only CSS is currently hot reloaded
    new webpack.HotModuleReplacementPlugin(),
    new webpack.NoEmitOnErrorsPlugin()
  ]
};
