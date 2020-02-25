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

beforeEach(async (done) => {
  await page.goto(appUrlBase);
  done();
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

  beforeEach(async (done) => {
    if (config.withSubset) await subset();
    done(); // https://github.com/facebook/jest/issues/1256
  });

  test("create a category", async () => {
    await assertCategoryDoesNotExist("test-category");
    await cxgActions.createCategory("test-category");
    await assertCategoryExists("test-category");
  });

  test("delete a category", async () => {
    const categoryName = "cluster-test";
    await assertCategoryExists(categoryName);
    await cxgActions.deleteCategory(categoryName);
    await assertCategoryDoesNotExist(categoryName);
  });


  test("rename a category", async () => {
    const oldCategoryName = "cluster-test";
    const newCategoryName = "cluster-for-real";
    await assertCategoryExists(oldCategoryName);
    await cxgActions.renameCategory(oldCategoryName, newCategoryName);
    await assertCategoryDoesNotExist(oldCategoryName);
    await assertCategoryExists(newCategoryName);
  });

  test("create a label", async () => {
    const categoryName = "cluster-test";
    const labelName = "new-label";
    await assertLabelDoesNotExist(categoryName, labelName);
    await cxgActions.createLabel(categoryName, labelName);
    await assertLabelExists(categoryName, labelName);
  });

  test("delete a label", async () => {
    const categoryName = "cluster-test";
    const labelName = "three";
    await assertLabelExists(categoryName, labelName);
    await cxgActions.deleteLabel(categoryName, labelName);
    await assertLabelDoesNotExist(categoryName, labelName);
  });

  test("rename a label", async () => {
    const categoryName = "cluster-test";
    const oldLabelName = "two";
    const newLabelName = "bazillion";
    await assertLabelExists(categoryName, oldLabelName);
    await cxgActions.renameLabel(categoryName, oldLabelName, newLabelName);
    await assertLabelDoesNotExist(categoryName, oldLabelName);
    await assertLabelExists(categoryName, newLabelName);
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

  test("undo/redo category creation", async () => {
    const categoryName = "heisenbug-category";
    await assertCategoryDoesNotExist(categoryName);
    await cxgActions.createCategory(categoryName);
    await assertCategoryExists(categoryName);
    await utils.clickOn("undo");
    await assertCategoryDoesNotExist(categoryName);
    await utils.clickOn("redo");
    await assertCategoryExists(categoryName);
  });

  test("undo/redo category deletion", async () => {
    const categoryName = "cluster-test";
    await assertCategoryExists(categoryName);
    await cxgActions.deleteCategory(categoryName);
    await assertCategoryDoesNotExist(categoryName);
    await utils.clickOn("undo");
    await assertCategoryExists(categoryName);
    await utils.clickOn("redo");
    await assertCategoryDoesNotExist(categoryName);
  });

  test("undo/redo category rename", async () => {
    const oldCategoryName = "cluster-test";
    const newCategoryName = "heisenbug-category";
    await assertCategoryExists(oldCategoryName);
    await assertCategoryDoesNotExist(newCategoryName);
    await cxgActions.renameCategory(oldCategoryName, newCategoryName);
    await assertCategoryExists(newCategoryName);
    await assertCategoryDoesNotExist(oldCategoryName);
    await utils.clickOn("undo");
    await assertCategoryExists(oldCategoryName);
    await assertCategoryDoesNotExist(newCategoryName);
    await utils.clickOn("redo");
    await assertCategoryExists(newCategoryName);
    await assertCategoryDoesNotExist(oldCategoryName);
  });

  test("undo/redo label creation", async () => {
    const categoryName = "cluster-test";
    const labelName = "heisenbug-label";
    await assertLabelDoesNotExist(categoryName, labelName);
    await cxgActions.createLabel(categoryName, labelName);
    await assertLabelExists(categoryName, labelName);
    await utils.clickOn("undo");
    await assertLabelDoesNotExist(categoryName);
    await utils.clickOn("redo");
    await assertLabelExists(categoryName, labelName);
  });

  test("undo/redo label deletion", async () => {
    const categoryName = "cluster-test";
    const labelName = "five";
    await assertLabelExists(categoryName, labelName);
    await cxgActions.deleteLabel(categoryName, labelName);
    await assertLabelDoesNotExist(categoryName);
    await utils.clickOn("undo");
    await assertLabelExists(categoryName, labelName);
    await utils.clickOn("redo");
    await assertLabelDoesNotExist(categoryName);
  });

  test("undo/redo label rename", async () => {
    const categoryName = "cluster-test";
    const oldLabelName = "six";
    const newLabelName = "schroedinger-label";
    await assertLabelExists(categoryName, oldLabelName);
    await assertLabelDoesNotExist(categoryName, newLabelName);
    await cxgActions.renameLabel(categoryName, oldLabelName, newLabelName);
    await assertLabelExists(categoryName, newLabelName);
    await assertLabelDoesNotExist(categoryName, oldLabelName);
    await utils.clickOn("undo");
    await assertLabelExists(categoryName, oldLabelName);
    await assertLabelDoesNotExist(categoryName, newLabelName);
    await utils.clickOn("redo");
    await assertLabelExists(categoryName, newLabelName);
    await assertLabelDoesNotExist(categoryName, oldLabelName);
  });

  async function assertCategoryExists(categoryName) {
    const result = await utils.waitByID(`${categoryName}:category-expand`);
    expect(await result.evaluate(node => node.innerText)).toBe(categoryName);
  }

  async function assertCategoryDoesNotExist(categoryName) {
    const result = await page.$(`[data-testid='${categoryName}:category-expand']`);
    expect(result).toBeNull();
  }

  async function assertLabelExists(categoryName, labelName) {
    const category = await utils.waitByID(`${categoryName}:category-expand`);
    expect(category).not.toBeNull();
    await cxgActions.expandCategory(categoryName);
    const previous = await utils.waitByID(`categorical-value-${categoryName}-${labelName}`);
    expect(await previous.evaluate(node => node.innerText)).toBe(labelName);
  }

  async function assertLabelDoesNotExist(categoryName, labelName) {
    await cxgActions.expandCategory(categoryName);
    const result = await page.$(`[data-testid='categorical-value-${categoryName}-${labelName}']`);
    expect(result).toBeNull();
  }
});
