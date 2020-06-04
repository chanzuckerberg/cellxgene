/*
Tests included in this file are specific to annotation features
*/
import { appUrlBase, DATASET } from "./config";
import { setupTestBrowser } from "./testBrowser";
import { datasets } from "./data";

let browser;
let page;
let utils;
let actions;
const data = datasets[DATASET];

beforeAll(async () => {
  [browser, page, utils, actions] = await setupTestBrowser();
});

afterAll(() => {
  if (browser !== undefined) browser.close();
});

describe.each([
  { withSubset: true, tag: "subset" },
  { withSubset: false, tag: "whole" },
])("annotations", (config) => {
  const perTestCategoryName = "per-test-category";
  const perTestLabelName = "per-test-label";

  beforeEach(async () => {
    await page.goto(appUrlBase);
    // wait for the page to load
    await utils.waitByClass("autosave-complete");
    // setup the test fixtures
    await actions.createCategory(perTestCategoryName);
    await actions.createLabel(perTestCategoryName, perTestLabelName);
    if (config.withSubset)
      await actions.subset({ x1: 0.1, y1: 0.1, x2: 0.8, y2: 0.8 });
    await utils.waitByClass("autosave-complete");
  });

  afterEach(async () => {
    await deleteCategoryIfExists(perTestCategoryName);
    await utils.waitByClass("autosave-complete");
  });

  test("create a category", async () => {
    const categoryName = `category-created-${config.tag}`;
    await assertCategoryDoesNotExist(categoryName);
    await actions.createCategory(categoryName);
    await assertCategoryExists(categoryName);
  });

  test("delete a category", async () => {
    await actions.deleteCategory(perTestCategoryName);
    await assertCategoryDoesNotExist(perTestCategoryName);
  });

  test("rename a category", async () => {
    const newCategoryName = `cluster-for-real-${config.tag}`;
    await actions.renameCategory(perTestCategoryName, newCategoryName);
    await assertCategoryDoesNotExist(perTestCategoryName);
    await assertCategoryExists(newCategoryName);
  });

  test("create a label", async () => {
    const labelName = `new-label-${config.tag}`;
    await assertLabelDoesNotExist(perTestCategoryName, labelName);
    await actions.createLabel(perTestCategoryName, labelName);
    await assertLabelExists(perTestCategoryName, labelName);
  });

  test("delete a label", async () => {
    await actions.deleteLabel(perTestCategoryName, perTestLabelName);
    await assertLabelDoesNotExist(perTestCategoryName, perTestLabelName);
  });

  test("rename a label", async () => {
    const newLabelName = "my-cool-new-label";
    await assertLabelDoesNotExist(perTestCategoryName, newLabelName);
    await actions.renameLabel(
      perTestCategoryName,
      perTestLabelName,
      newLabelName
    );
    await assertLabelDoesNotExist(perTestCategoryName, perTestLabelName);
    await assertLabelExists(perTestCategoryName, newLabelName);
  });

  test("check cell count for a label loaded from file", async () => {
    const categoryName = "cluster-test";
    const labelName = "four";
    await actions.expandCategory(categoryName);
    const result = await utils.waitByID(
      `categorical-value-count-${categoryName}-${labelName}`
    );
    expect(await result.evaluate((node) => node.innerText)).toBe(
      data.annotationsFromFile.count.bySubsetConfig[config.withSubset]
    );
  });

  test("assign cells to a label", async () => {
    await actions.expandCategory(perTestCategoryName);
    const lassoSelection = await actions.calcDragCoordinates(
      "layout-graph",
      data.categoryLabel.lasso["coordinates-as-percent"]
    );
    await actions.drag(
      "layout-graph",
      lassoSelection.start,
      lassoSelection.end,
      true
    );
    await utils.waitByID("lasso-element", { visible: true });
    await utils.clickOn(
      `${perTestCategoryName}:${perTestLabelName}:see-actions`
    );
    await utils.clickOn(
      `${perTestCategoryName}:${perTestLabelName}:add-current-selection-to-this-label`
    );
    const result = await utils.waitByID(
      `categorical-value-count-${perTestCategoryName}-${perTestLabelName}`
    );
    expect(await result.evaluate((node) => node.innerText)).toBe(
      data.categoryLabel.newCount.bySubsetConfig[config.withSubset]
    );
  });

  test("undo/redo category creation", async () => {
    const categoryName = `category-created-undo-${config.tag}`;
    await assertCategoryDoesNotExist(categoryName);
    await actions.createCategory(categoryName);
    await assertCategoryExists(categoryName);
    await utils.clickOn("undo");
    await assertCategoryDoesNotExist(categoryName);
    await utils.clickOn("redo");
    await assertCategoryExists(categoryName);
  });

  test("undo/redo category deletion", async () => {
    const categoryName = `category-deleted-undo-${config.tag}`;
    await actions.createCategory(categoryName);
    await assertCategoryExists(categoryName);
    await actions.deleteCategory(categoryName);
    await assertCategoryDoesNotExist(categoryName);
    await utils.clickOn("undo");
    await assertCategoryExists(categoryName);
    await utils.clickOn("redo");
    await assertCategoryDoesNotExist(categoryName);
  });

  test("undo/redo category rename", async () => {
    const newCategoryName = `category-renamed-undo-${config.tag}`;
    await assertCategoryDoesNotExist(newCategoryName);
    await actions.renameCategory(perTestCategoryName, newCategoryName);
    await assertCategoryExists(newCategoryName);
    await assertCategoryDoesNotExist(perTestCategoryName);
    await utils.clickOn("undo");
    await assertCategoryExists(perTestCategoryName);
    await assertCategoryDoesNotExist(newCategoryName);
    await utils.clickOn("redo");
    await assertCategoryExists(newCategoryName);
    await assertCategoryDoesNotExist(perTestCategoryName);
  });

  test("undo/redo label creation", async () => {
    const labelName = `label-created-undo-${config.tag}`;
    await assertLabelDoesNotExist(perTestCategoryName, labelName);
    await actions.createLabel(perTestCategoryName, labelName);
    await assertLabelExists(perTestCategoryName, labelName);
    await utils.clickOn("undo");
    await assertLabelDoesNotExist(perTestCategoryName);
    await utils.clickOn("redo");
    await assertLabelExists(perTestCategoryName, labelName);
  });

  test("undo/redo label deletion", async () => {
    await actions.deleteLabel(perTestCategoryName, perTestLabelName);
    await assertLabelDoesNotExist(perTestCategoryName);
    await utils.clickOn("undo");
    await assertLabelExists(perTestCategoryName, perTestLabelName);
    await utils.clickOn("redo");
    await assertLabelDoesNotExist(perTestCategoryName);
  });

  test("undo/redo label rename", async () => {
    const newLabelName = `label-renamed-undo-${config.tag}`;
    await assertLabelDoesNotExist(perTestCategoryName, newLabelName);
    await actions.renameLabel(
      perTestCategoryName,
      perTestLabelName,
      newLabelName
    );
    await assertLabelExists(perTestCategoryName, newLabelName);
    await assertLabelDoesNotExist(perTestCategoryName, perTestLabelName);
    await utils.clickOn("undo");
    await assertLabelExists(perTestCategoryName, perTestLabelName);
    await assertLabelDoesNotExist(perTestCategoryName, newLabelName);
    await utils.clickOn("redo");
    await assertLabelExists(perTestCategoryName, newLabelName);
    await assertLabelDoesNotExist(perTestCategoryName, perTestLabelName);
  });

  async function assertCategoryExists(categoryName) {
    const handle = await utils.waitByID(`${categoryName}:category-label`);

    const result = await handle.evaluate((node) =>
      node.getAttribute("aria-label")
    );

    expect(result).toBe(categoryName);
  }

  async function assertCategoryDoesNotExist(categoryName) {
    const result = await page.$(
      `[data-testid='${categoryName}:category-label']`
    );
    expect(result).toBeNull();
  }

  async function assertLabelExists(categoryName, labelName) {
    const category = await utils.waitByID(`${categoryName}:category-expand`);
    expect(category).not.toBeNull();
    await actions.expandCategory(categoryName);
    const previous = await utils.waitByID(
      `categorical-value-${categoryName}-${labelName}`
    );
    expect(
      await previous.evaluate((node) => node.getAttribute("aria-label"))
    ).toBe(labelName);
  }

  async function assertLabelDoesNotExist(categoryName, labelName) {
    await actions.expandCategory(categoryName);
    const result = await page.$(
      `[data-testid='categorical-value-${categoryName}-${labelName}']`
    );
    expect(result).toBeNull();
  }

  async function deleteCategoryIfExists(categoryName) {
    try {
      const category = await page.waitForSelector(
        `[data-testid='${categoryName}:category-expand']`,
        { timeout: 200 }
      );
      if (category !== null) return await actions.deleteCategory(categoryName);
    } catch {}
    return null;
  }
});
