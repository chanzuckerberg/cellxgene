module.exports = {
  babelrc: false,
  presets: [
    [
      "@babel/preset-env",
      {
        useBuiltIns: "entry",
        corejs: 3,
        modules: false,
      },
    ],
    "@babel/preset-react",
  ],
  plugins: [
    "@babel/plugin-proposal-function-bind",
    ["@babel/plugin-proposal-decorators", { legacy: true }],
    ["@babel/plugin-proposal-class-properties", { loose: true }],
    ["@babel/plugin-transform-private-methods", { loose: true }],
    ["@babel/plugin-transform-private-property-in-object", { loose: true }],
    "@babel/plugin-proposal-export-namespace-from",
    "@babel/plugin-transform-react-constant-elements",
    "@babel/plugin-transform-runtime",
    "@babel/plugin-proposal-optional-chaining",
    "@babel/plugin-proposal-nullish-coalescing-operator",
  ],
};
