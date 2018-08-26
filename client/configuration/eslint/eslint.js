module.exports = {
  root: true,
  parser: "babel-eslint",
  extends: "airbnb",
  env: { browser: true, commonjs: true, es6: true },
  globals: { expect: true },
  parserOptions: {
    ecmaVersion: 2017,
    sourceType: "module",
    ecmaFeatures: {
      jsx: true,
      generators: true
    }
  },
  rules: {
    "no-magic-numbers": "off",
    "func-style": "off",
    "arrow-parens": "off",
    "no-use-before-define": "off",
    "react/jsx-filename-extension": "off",
    "comma-dangle": "off",
    "no-underscore-dangle": "off",
    quotes: ["error", "double"],
    "implicit-arrow-linebreak": "off",
    "operator-linebreak": ["error", "after"]
  }
};
