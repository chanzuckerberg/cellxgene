module.exports = {
  root: true,
  parser: "babel-eslint",
  extends: [
    "airbnb",
    "plugin:eslint-comments/recommended",
    "plugin:@blueprintjs/recommended",
    "plugin:compat/recommended",
    "plugin:prettier/recommended",
    "prettier/react",
  ],
  settings: {
    polyfills: [
      "TextDecoder",
      "TextEncoder",
      "fetch",
      "Request",
      "Response",
      "Headers",
      "AbortController",
    ],
  },
  env: { browser: true, commonjs: true, es6: true },
  globals: {
    expect: true,
    jest: true,
    jestPuppeteer: true,
    it: true,
    page: true,
    browser: true,
    context: true,
    beforeEach: true,
  },
  parserOptions: {
    ecmaVersion: 2017,
    sourceType: "module",
    ecmaFeatures: {
      jsx: true,
      generators: true,
    },
  },
  rules: {
    "react/jsx-no-target-blank": "off",
    "eslint-comments/require-description": ["error"],
    "no-magic-numbers": "off",
    "no-nested-ternary": "off",
    "func-style": "off",
    "arrow-parens": "off",
    "no-use-before-define": "off",
    "react/jsx-filename-extension": "off",
    "comma-dangle": "off",
    "no-underscore-dangle": "off",
    "implicit-arrow-linebreak": "off",
    "no-console": "off",
    "spaced-comment": ["error", "always", { exceptions: ["*"] }],
    "no-param-reassign": "off",
    "object-curly-newline": ["error", { consistent: true }],
    "react/prop-types": [0],
    "space-before-function-paren": "off",
    "function-paren-newline": "off",
    "prefer-destructuring": ["error", { object: true, array: false }],
    "import/prefer-default-export": "off",
    "no-restricted-syntax": [
      "error",
      "ForInStatement",
      "LabeledStatement",
      "WithStatement",
    ],
    "import/no-extraneous-dependencies": [
      "error",
      {
        devDependencies: true,
      },
    ],
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
