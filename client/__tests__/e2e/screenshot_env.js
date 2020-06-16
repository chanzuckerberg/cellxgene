const PuppeteerEnvironment = require("jest-environment-puppeteer");
require("jest-circus");

const takeScreenshot = require("./takeScreenshot");

// (thuang): Please make sure this number matches
// `jest.retryTimes(N)` in `puppeteer.setup.js`
const RETRY_ATTEMPTS = 2;

class ScreenshotEnvironment extends PuppeteerEnvironment {
  async handleTestEvent(event, state) {
    if (event.name === "error") {
      console.log("error", JSON.stringify(event));
    }

    if (event.name === "test_fn_failure" || event.name === "hook_failure") {
      // (thuang): We only want to take screenshot on the last try
      if (state.currentlyRunningTest.invocations <= RETRY_ATTEMPTS) return;

      console.log("===> Failure event\n", new Date(), event);

      await takeScreenshot(state.currentlyRunningTest.name, this.global.page);
    }
  }
}

module.exports = ScreenshotEnvironment;
