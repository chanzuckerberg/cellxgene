import puppeteer from "puppeteer";
import { DEBUG, DEV } from "./config";
import { puppeteerUtils } from "./puppeteerUtils";
import { cellxgeneActions } from "./cellxgeneActions";

export async function setupTestBrowser() {
  const browserViewport = { width: 1280, height: 960 };
  const browserParams = DEV
    ? {
        headless: false,
        slowMo: 5,
        args: [
          `--window-size=${browserViewport.width},${browserViewport.height}`,
        ],
      }
    : DEBUG
    ? {
        headless: false,
        slowMo: 100,
        devtools: true,
        args: [
          `--window-size=${browserViewport.width + 560},${
            browserViewport.height
          }`,
        ],
      }
    : {
        args: [
          `--window-size=${browserViewport.width},${browserViewport.height}`,
        ],
      };
  const browser = await puppeteer.launch(browserParams);
  const page = await browser.pages().then((pages) => pages[0]);
  await page.setViewport(browserViewport);
  if (DEV || DEBUG) {
    page.on("console", async (msg) => {
      // If there is a console.error but an error is not thrown, this will ensure the test fails
      if (msg.type() === "error") {
        const errorMsgText = await Promise.all(
          // TODO can we do this without internal properties?
          msg.args().map((arg) => arg._remoteObject.description)
        );
        throw new Error(`Console error: ${errorMsgText}`);
      }
      console.log(`PAGE LOG: ${msg.text()}`);
    });
  }
  page.on("pageerror", (err) => {
    throw new Error(`Console error: ${err}`);
  });
  const utils = puppeteerUtils(page);
  const cxgActions = cellxgeneActions(page, utils);
  return [browser, page, utils, cxgActions];
}
