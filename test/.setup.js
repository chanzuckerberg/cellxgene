var babelConfig = require('../configuration/babel/babel.test.js');
require('babel-register')(babelConfig);

// Deobfuscate CSS modules classes for testing
require('css-modules-require-hook')({ generateScopedName: '[local]' });

require.extensions['.svg'] = () => null;

var jsdom = require('jsdom').jsdom;

var chai = require('chai');
var sinonChai = require('sinon-chai');

var exposedProperties = [ 'window', 'navigator', 'document' ];

global.document = jsdom('');
global.window = document.defaultView;

global.expect = chai.expect;

chai.use(sinonChai);

Object.keys(document.defaultView).forEach(property => {
  if (typeof global[property] === 'undefined') {
    exposedProperties.push(property);
    global[property] = document.defaultView[property];
  }
});

global.navigator = { userAgent: 'node.js' };

documentRef = document;
