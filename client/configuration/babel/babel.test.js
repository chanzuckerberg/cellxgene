module.exports = {
  babelrc: false,
  presets: [
    ["modern-browsers", { loose: true, modules: false }],
    "stage-0",
    "react"
  ],
  plugins: [
    "babel-plugin-transform-react-constant-elements",
    "babel-plugin-transform-decorators-legacy",
    "babel-plugin-transform-runtime"
  ]
};
