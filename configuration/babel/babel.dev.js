module.exports = {
  babelrc: false,
  cacheDirectory: true,
  presets: [
    [ "env", { loose: true, modules: false } ],
    "stage-0",
    "react"
  ],
  plugins: [
    "react-hot-loader/babel",
    "babel-plugin-transform-decorators-legacy"
  ]
};
