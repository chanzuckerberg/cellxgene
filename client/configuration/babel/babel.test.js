module.exports = {
  babelrc: false,
  presets: [
    ["modern-browsers", { loose: true, modules: false }],
    "stage-0",
    "react"
  ],
  plugins: ["istanbul"]
};
