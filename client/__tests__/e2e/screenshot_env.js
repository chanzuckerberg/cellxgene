const PuppeteerEnvironment = require("jest-environment-puppeteer");
require("jest-circus");
const ENV_DEFAULT = require("../../../environment.default.json");

const takeScreenshot = require("./takeScreenshot");

class ScreenshotEnvironment extends PuppeteerEnvironment {
  async handleTestEvent(event, state) {
    if (["test_start", "test_done"].includes(event.name)) {
      console.log("------------------event name:\n", event.name);
      console.log("~~~~ Current test errors\n", new Date(), event.test.errors);
      console.log("~~~~ Current test\n", new Date(), event.test);
    }

    if (event.name === "error") {
      console.log("error event:", JSON.stringify(event));
    }

    if (event.name === "test_fn_failure" || event.name === "hook_failure") {
      console.log("------------------event name:\n", event.name);
      console.log(">>>> Current state\n", new Date(), state);
      console.log("===> Failure event\n", new Date(), event);

      // (thuang): We only want to take screenshot on the last try
      if (
        state.currentlyRunningTest.invocations <= ENV_DEFAULT.RETRY_ATTEMPTS
      ) {
        return;
      }

      await takeScreenshot(state.currentlyRunningTest.name, this.global.page);
    }
  }
}

module.exports = ScreenshotEnvironment;
