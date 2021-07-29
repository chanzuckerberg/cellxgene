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
    // eslint-disable-next-line prettier/prettier --- FIXME: disabled temporarily on migrate to TS.
    "@babel/preset-typescript"
  ],
  plugins: [
    "@babel/plugin-proposal-function-bind",
    ["@babel/plugin-proposal-decorators", { legacy: true }],
    ["@babel/plugin-proposal-class-properties", { loose: true }],
    "@babel/plugin-proposal-export-namespace-from",
    "@babel/plugin-transform-react-constant-elements",
    "@babel/plugin-transform-runtime",
    "@babel/plugin-proposal-optional-chaining",
    "@babel/plugin-proposal-nullish-coalescing-operator",
  ],
};
