module.exports = {
  babelrc: false,
  presets: [
    [ 'es2015', { loose: true, modules: false } ],
    'stage-0',
    'react'
  ],
  plugins: [ 'babel-plugin-transform-react-constant-elements' ]
    .map(require.resolve)
    .concat([ [ require.resolve('babel-plugin-transform-runtime') ] ])
};
