/**
 * `client/jest-puppeteer.config.js` is for configuring Puppeteer's launch config options
 * `client/__tests__/e2e/puppeteer.setup.js` is for configuring `jest`, `browser`,
 * and `page` objects
 */

const ENV_DEFAULT = require("../environment.default.json");

const jestEnv = process.env.JEST_ENV || ENV_DEFAULT.JEST_ENV;
// const isHeadful =
//   process.env.HEADFUL === "true" || process.env.HEADLESS === "false";

const DEFAULT_LAUNCH_CONFIG = {
  headless: true,
  args: ["--ignore-certificate-errors", "--ignore-ssl-errors"],
  ignoreHTTPSErrors: true,
  defaultViewport: {
    width: 1280,
    height: 960,
  },
};

const LAUNCH_CONFIG_BY_ENV = {
  [ENV_DEFAULT.DEBUG]: {
    ...DEFAULT_LAUNCH_CONFIG,
    headless: true,
    slowMo: 100,
    devtools: true,
    defaultViewport: {
      width: DEFAULT_LAUNCH_CONFIG.defaultViewport.width,
      height: DEFAULT_LAUNCH_CONFIG.defaultViewport.height + 560,
    },
  },
  [ENV_DEFAULT.DEV]: {
    ...DEFAULT_LAUNCH_CONFIG,
    headless: true,
    slowMo: 5,
  },
};

const launchConfig = LAUNCH_CONFIG_BY_ENV[jestEnv] || DEFAULT_LAUNCH_CONFIG;

module.exports = {
  browserContext: "incognito",
  launch: launchConfig,
};
