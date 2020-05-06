/*
NOT run in Travis CI

UX tests using puppeteer to be run locally.

To run locally, ensure you are running the client is running on port 3000.
Then run jest --verbose false --config __tests__/e2e/e2eJestConfig.json feature.
 */

import puppeteer from "puppeteer";
import { appUrlBase, DEBUG, DEV, DATASET } from "./config";
import { puppeteerUtils, cellxgeneActions } from "./puppeteerUtils";
import { datasets } from "./data";

let browser;
let page;
let utils;
let cxgActions;
let spy;
const browserViewport = { width: 1280, height: 960 };
const data = datasets[DATASET].features;

if (DEBUG) jest.setTimeout(100000);
if (DEV) jest.setTimeout(10000);

beforeAll(async () => {
  const browserParams = DEV
    ? { headless: false, slowMo: 5 }
    : DEBUG
    ? { headless: false, slowMo: 100, devtools: true }
    : {};
  browser = await puppeteer.launch(browserParams);
  page = await browser.newPage();
  await page.setViewport(browserViewport);
  if (DEV || DEBUG) {
    page.on("console", (msg) => console.log(`PAGE LOG: ${msg.text()}`));
  }
  page.on("pageerror", (err) => {
    throw new Error(`Console error: ${err}`);
  });
  utils = puppeteerUtils(page);
  cxgActions = cellxgeneActions(page);
});

beforeEach(async () => {
  await page.goto(appUrlBase);
});

afterAll(() => {
  if (!DEBUG) {
    browser.close();
  }
});

describe("zoom interaction", async () => {
  // Skip this test since UI is to hide lasso path when switching modes
  test.skip("lasso visible after switching modes to pan/zoom", async () => {
    const lassoSelection = await cxgActions.calcDragCoordinates(
      "layout-graph",
      data.panzoom.lasso["coordinates-as-percent"]
    );
    await cxgActions.drag(
      "layout-graph",
      lassoSelection.start,
      lassoSelection.end,
      true
    );
    await utils.waitByID("lasso-element", { visible: true });
    await utils.clickOn("mode-pan-zoom");
    await utils.waitByID("lasso-element", { visible: true });
  });

  test("pan zoom mode resets lasso selection", async () => {
    const lassoSelection = await cxgActions.calcDragCoordinates(
      "layout-graph",
      data.panzoom.lasso["coordinates-as-percent"]
    );
    await cxgActions.drag(
      "layout-graph",
      lassoSelection.start,
      lassoSelection.end,
      true
    );
    await utils.waitByID("lasso-element", { visible: true });
    const initialCount = await cxgActions.cellSet(1);
    expect(initialCount).toBe(data.panzoom.lasso.count);
    await utils.clickOn("mode-pan-zoom");
    await utils.clickOn("mode-lasso");
    const modeSwitchCount = await cxgActions.cellSet(1);
    expect(modeSwitchCount).toBe(initialCount);
  });

  test("lasso moves after pan", async () => {
    const lassoSelection = await cxgActions.calcDragCoordinates(
      "layout-graph",
      data.panzoom.lasso["coordinates-as-percent"]
    );
    await cxgActions.drag(
      "layout-graph",
      lassoSelection.start,
      lassoSelection.end,
      true
    );
    await utils.waitByID("lasso-element", { visible: true });
    const initialCount = await cxgActions.cellSet(1);
    expect(initialCount).toBe(data.panzoom.lasso.count);
    await utils.clickOn("mode-pan-zoom");
    const panCoords = await cxgActions.calcDragCoordinates(
      "layout-graph",
      data.panzoom.lasso["coordinates-as-percent"]
    );
    await cxgActions.drag(
      "layout-graph",
      panCoords.start,
      panCoords.end,
      false
    );
    await utils.clickOn("mode-lasso");
    const panCount = await cxgActions.cellSet(2);
    expect(panCount).toBe(initialCount);
  });
});
