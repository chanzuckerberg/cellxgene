/*
Tests included in this file are specific to annotation features
*/
import {appUrlBase, DEBUG, DEV, DATASET} from "./config";
import {setupTestBrowser} from "./puppeteerUtils";
import {datasets} from "./data";

let browser, page, utils, cxgActions;
const browserViewport = {width: 1280, height: 960};
const data = datasets[DATASET];

beforeAll(async () => {
  [browser, page, utils, cxgActions] = await setupTestBrowser(browserViewport);
});

beforeEach(async () => {
  await page.goto(appUrlBase);
});

afterAll(() => {
  if (!DEBUG) browser.close();
});

describe("did launch", () => {
  test("page launched", async () => {
    let el = await utils.getOneElementInnerHTML("[data-testid='header']");
    expect(el).toBe(data.title);
  });
});


describe.each([
  {withSubset: true},
  {withSubset: false}
])("annotations", (config) => {

  async function subset() {
    const lassoSelection = await cxgActions.calcDragCoordinates(
      "layout-graph",
      { x1: 0.10, y1: 0.10, x2: 0.80, y2: 0.80 }
    );
    await cxgActions.drag(
      "layout-graph",
      lassoSelection.start,
      lassoSelection.end,
      true
    );
    await utils.clickOn("subset-button");
    const coordinate = await cxgActions.calcCoordinate("layout-graph", 0.9, 0.9);
    await cxgActions.clickOnCoordinate("layout-graph", coordinate);
  }

  beforeEach(async () => {
    if (config.withSubset) await subset();
  });

  test("create a category", async () => {
    await utils.clickOn("open-annotation-dialog");
    await utils.typeInto("new-category-name", "test-category-name");
    await utils.clickOn("submit-category");
    const result = await utils.waitByID("test-category-name:category-expand");
    expect(await result.evaluate(node => node.innerText)).toBe("test-category-name");
  });

  test("delete a category", async () => {
    const previous = await utils.waitByID("cluster-test:category-expand");
    expect(await previous.evaluate(node => node.innerText)).toBe("cluster-test");
    await utils.hoverOn("cluster-test:see-actions");
    await utils.clickOn("cluster-test:delete-category");
    const result = await page.$("[data-testid='cluster-test:category-expand']");
    expect(result).toBeNull();
  });

  test("rename a category", async () => {
    const previous = await utils.waitByID("cluster-test:category-expand");
    expect(await previous.evaluate(node => node.innerText)).toBe("cluster-test");
    await utils.hoverOn("cluster-test:see-actions");
    await utils.clickOn("cluster-test:edit-category-mode");
    await utils.typeInto("cluster-test:edit-category-name-text", "-renamed");
    await utils.clickOn("cluster-test:submit-category-edit");
    const result = await utils.waitByID("cluster-test-renamed:category-expand");
    expect(await result.evaluate(node => node.innerText)).toBe("cluster-test-renamed");
  });

  test("create a label", async () => {
    await utils.hoverOn("cluster-test:see-actions");
    await utils.clickOn("cluster-test:add-new-label-to-category");
    await utils.typeInto("cluster-test:new-label-name", "test-label-name");
    await utils.clickOn("cluster-test:submit-label");
    await cxgActions.expandCategory("cluster-test");
    const result = await utils.waitByID("categorical-value-cluster-test-test-label-name");
    expect(await result.evaluate(node => node.innerText)).toBe("test-label-name");
  });

  test("delete a label", async () => {
    await cxgActions.expandCategory("cluster-test");
    const previous = await utils.waitByID("categorical-value-cluster-test-three");
    expect(await previous.evaluate(node => node.innerText)).toBe("three");
    await utils.hoverOn("cluster-test:three:see-actions");
    await utils.clickOn("cluster-test:three:delete-label");
    const result = await page.$("[data-testid='categorical-value-cluster-test-three']");
    expect(result).toBeNull();
  });

  test("rename a label", async () => {
    await cxgActions.expandCategory("cluster-test");
    const previous = await utils.waitByID("categorical-value-cluster-test-four");
    expect(await previous.evaluate(node => node.innerText)).toBe("four");
    await utils.hoverOn("cluster-test:four:see-actions");
    await utils.clickOn("cluster-test:four:edit-label");
    await utils.typeInto("cluster-test:four:edit-label-name", ".");
    await utils.clickOn("cluster-test:four:submit-label-edit");
    const result = await utils.waitByID("categorical-value-cluster-test-four.");
    expect(await result.evaluate(node => node.innerText)).toBe("four.");
  });

  test("assign cells to a label", async () => {
    await cxgActions.expandCategory("cluster-test");

    const lassoSelection = await cxgActions.calcDragCoordinates(
      "layout-graph",
      data.categoryLabel.lasso["coordinates-as-percent"]
    );
    await cxgActions.drag(
      "layout-graph",
      lassoSelection.start,
      lassoSelection.end,
      true
    );
    await utils.waitByID("lasso-element", {visible: true});
    const initialCount = await cxgActions.cellSet(1);
    expect(initialCount).toBe(data.categoryLabel.lasso.count);
    await utils.hoverOn("cluster-test:one:see-actions");
    await utils.clickOn("cluster-test:one:add-current-selection-to-this-label");
    const result = await utils.waitByID("categorical-value-count-cluster-test-one");
    expect(await result.evaluate(node => node.innerText)).toBe(
      data.categoryLabel.newCount.bySubsetConfig[config.withSubset]
    );
  });
});
