module.exports = {
  root: true,
  parser: 'babel-eslint',
  extends: 'formidable/configurations/es6-react',
  env: { browser: true, commonjs: true, es6: true, node: true, mocha: true },
  globals: { expect: true },
  parserOptions: {
    ecmaVersion: 6,
    sourceType: 'module',
    ecmaFeatures: {
      jsx: true,
      generators: true,
      experimentalObjectRestSpread: true
    }
  },
  rules: {
    quotes: [ 2, 'single', { allowTemplateLiterals: true } ],
    'no-magic-numbers': 'off',
    'func-style': 'off',
    'arrow-parens': 'off',
    'no-use-before-define': 'off',
    'react/jsx-filename-extension': 'off',
    'react/require-extension': 'off',
    'react/no-multi-comp': 'warn',
    'react/prop-types': 'warn',
    'react/sort-comp': 'warn',
    'react/sort-prop-types': 'warn'
  }
};
