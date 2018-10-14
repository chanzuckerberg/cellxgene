module.exports = {
  babelrc: false,
  presets: [
    ["modern-browsers", { loose: true, modules: false }],
    "@babel/preset-react"
  ],
  plugins: [
    "@babel/plugin-proposal-function-bind",
    "@babel/plugin-proposal-class-properties",
    ["@babel/plugin-proposal-decorators", { legacy: true }],
    "@babel/plugin-proposal-export-namespace-from",
    "@babel/plugin-transform-react-constant-elements",
    "@babel/plugin-transform-runtime"
  ]
};
