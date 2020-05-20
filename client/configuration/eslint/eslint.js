module.exports = {
  root: true,
  parser: "babel-eslint",
  extends: ["airbnb", "plugin:prettier/recommended", "prettier/react"],
  env: { browser: true, commonjs: true, es6: true },
  globals: { expect: true },
  parserOptions: {
    ecmaVersion: 2017,
    sourceType: "module",
    ecmaFeatures: {
      jsx: true,
      generators: true,
    },
  },
  rules: {
    "no-magic-numbers": "off",
    "no-nested-ternary": "off",
    "func-style": "off",
    "arrow-parens": "off",
    "no-use-before-define": "off",
    "react/jsx-filename-extension": "off",
    "comma-dangle": "off",
    "no-underscore-dangle": "off",
    quotes: ["error", "double"],
    "implicit-arrow-linebreak": "off",
    "operator-linebreak": [
      "error",
      "after",
      { overrides: { "?": "before", ":": "before" } },
    ],
    "no-console": "off",
    "spaced-comment": ["error", "always", { exceptions: ["*"] }],
    "no-param-reassign": "off",
    "object-curly-newline": ["error", { consistent: true }],
    "react/prop-types": [0],
    "space-before-function-paren": "off",
    "function-paren-newline": "off",
    "prefer-destructuring": ["error", { object: true, array: false }],
  },
  overrides: [
    {
      files: ["**/*.test.js"],
      env: {
        jest: true, // now **/*.test.js files' env has both es6 *and* jest
      },
      // Can't extend in overrides: https://github.com/eslint/eslint/issues/8813
      // "extends": ["plugin:jest/recommended"]
      plugins: ["jest"],
      rules: {
        "jest/no-disabled-tests": "warn",
        "jest/no-focused-tests": "error",
        "jest/no-identical-title": "error",
        "jest/prefer-to-have-length": "warn",
        "jest/valid-expect": "error",
      },
    },
  ],
};
