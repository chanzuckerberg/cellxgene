module.exports = {
  babelrc: false,
  cacheDirectory: true,
  presets: [
    [ 'es2015', { loose: true, modules: false } ],
    'stage-0',
    'react'
  ],
  plugins: [ 'react-hot-loader/babel' ]
};
