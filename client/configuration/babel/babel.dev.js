module.exports = {
  babelrc: false,
  cacheDirectory: true,
  presets: [
    ["modern-browsers", { loose: true, modules: false }],
    "@babel/preset-react"
  ],
  plugins: [
    "@babel/plugin-proposal-function-bind",
    ["@babel/plugin-proposal-decorators", { legacy: true }],
    ["@babel/plugin-proposal-class-properties", { loose: true }],
    "@babel/plugin-proposal-export-namespace-from",
    "@babel/plugin-proposal-optional-chaining",
    "@babel/plugin-proposal-nullish-coalescing-operator"
  ]
};
