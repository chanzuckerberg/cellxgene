module.exports = {
  babelrc: false,
  cacheDirectory: true,
  presets: [
    ["modern-browsers", { loose: true, modules: false }],
    "stage-0",
    "react"
  ],
  plugins: ["babel-plugin-transform-decorators-legacy"]
};
