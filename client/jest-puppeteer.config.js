/**
 * `client/jest-puppeteer.config.js` is for configuring Puppeteer's launch config options
 * `client/__tests__/e2e/puppeteer.setup.js` is for configuring `jest`, `browser`,
 * and `page` objects
 */

const DEBUG = "debug";
const DEV = "dev";
const PROD = "prod";
const jestEnv = process.env.JEST_ENV || PROD;
const isHeadful = process.env.HEADFUL === "true";

const DEFAULT_LAUNCH_CONFIG = {
  headless: !isHeadful,
  args: ["--ignore-certificate-errors", "--ignore-ssl-errors"],
  dumpio: true,
  ignoreHTTPSErrors: true,
  defaultViewport: {
    width: 1280,
    height: 960,
  },
};

const LAUNCH_CONFIG_BY_ENV = {
  [DEBUG]: {
    ...DEFAULT_LAUNCH_CONFIG,
    headless: false,
    slowMo: 100,
    devtools: true,
    defaultViewport: {
      width: DEFAULT_LAUNCH_CONFIG.defaultViewport.width,
      height: DEFAULT_LAUNCH_CONFIG.defaultViewport.height + 560,
    },
  },
  [DEV]: {
    ...DEFAULT_LAUNCH_CONFIG,
    headless: false,
    slowMo: 5,
  },
};

const launchConfig = LAUNCH_CONFIG_BY_ENV[jestEnv] || DEFAULT_LAUNCH_CONFIG;

module.exports = {
  browserContext: "incognito",
  launch: launchConfig,
};
