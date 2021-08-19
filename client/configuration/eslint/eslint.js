/* eslint-disable @blueprintjs/classes-constants -- we don't import blueprint here  */
module.exports = {
  root: true,
  parser: "@typescript-eslint/parser",
  extends: [
    "airbnb-typescript",
    "plugin:@typescript-eslint/recommended",
    "plugin:eslint-comments/recommended",
    "plugin:@blueprintjs/recommended",
    "plugin:compat/recommended",
    "plugin:prettier/recommended",
  ],
  settings: {
    // AbortController is not supported in iOS Safari 10.3, Chrome 61
    // Headers is not supported in iOS Safari 10.3
    polyfills: ["Headers", "AbortController"],
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
    project: "./tsconfig.json",
  },
  rules: {
    "react/jsx-no-target-blank": "off",
    "eslint-comments/require-description": ["error"],
    "no-magic-numbers": "off",
    "@typescript-eslint/no-magic-numbers": "off",
    "no-nested-ternary": "off",
    "func-style": "off",
    "arrow-parens": "off",
    "no-use-before-define": "off",
    "@typescript-eslint/no-use-before-define": "off",
    "react/jsx-filename-extension": "off",
    "comma-dangle": "off",
    "@typescript-eslint/comma-dangle": "off",
    "no-underscore-dangle": "off",
    // Override airbnb config to allow leading underscore
    // https://github.com/iamturns/eslint-config-airbnb-typescript/blob/master/lib/shared.js#L35
    "@typescript-eslint/naming-convention": [
      "error",
      {
        selector: "class",
        format: ["PascalCase"],
        leadingUnderscore: "allow",
      },
      {
        selector: "function",
        format: ["camelCase", "PascalCase"],
        leadingUnderscore: "allowSingleOrDouble",
      },
      {
        selector: "typeLike",
        format: ["PascalCase"],
      },
      {
        selector: "variable",
        format: ["camelCase", "PascalCase", "UPPER_CASE"],
        leadingUnderscore: "allowSingleOrDouble",
        trailingUnderscore: "allowDouble",
      },
    ],
    "implicit-arrow-linebreak": "off",
    "no-console": "off",
    "spaced-comment": ["error", "always", { exceptions: ["*"] }],
    "no-param-reassign": "off",
    "object-curly-newline": ["error", { consistent: true }],
    "react/prop-types": [0],
    "space-before-function-paren": "off",
    "@typescript-eslint/space-before-function-paren": "off",
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
      files: ["**/*.test.ts"],
      env: {
        jest: true, // now **/*.test.ts files' env has both es6 *and* jest
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
/* eslint-enable @blueprintjs/classes-constants -- we don't import blueprint here  */
