/*
Tests included in this file are specific to annotation features
*/
import { appUrlBase, DATASET } from "./config";
import { datasets } from "./data";

import {
  clickOn,
  goToPage,
  waitByClass,
  waitByID,
  getTestId,
  getTestClass,
  getAllByClass,
} from "./puppeteerUtils";

import {
  assertCategoryDoesNotExist,
  calcDragCoordinates,
  createCategory,
  createLabel,
  deleteCategory,
  deleteLabel,
  drag,
  expandCategory,
  renameCategory,
  renameLabel,
  subset,
  duplicateCategory,
} from "./cellxgeneActions";

const data = datasets[DATASET];

const perTestCategoryName = "TEST-CATEGORY";
const perTestLabelName = "TEST-LABEL";

async function setup(config) {
  await goToPage(appUrlBase);

  // setup the test fixtures
  await createCategory(perTestCategoryName);
  await createLabel(perTestCategoryName, perTestLabelName);

  if (config.withSubset) {
    await subset({ x1: 0.1, y1: 0.1, x2: 0.8, y2: 0.8 });
  }

  await waitByClass("autosave-complete");
}

describe.each([
  { withSubset: true, tag: "subset" },
  { withSubset: false, tag: "whole" },
])("annotations", (config) => {
  test("create a category", async () => {
    await setup(config);

    const categoryName = `category-created-${config.tag}`;

    await assertCategoryDoesNotExist(categoryName);

    await createCategory(categoryName);

    await assertCategoryExists(categoryName);
  });

  test("delete a category", async () => {
    await setup(config);

    await deleteCategory(perTestCategoryName);
    await assertCategoryDoesNotExist(perTestCategoryName);
  });

  test("rename a category", async () => {
    await setup(config);

    const newCategoryName = `NEW-${config.tag}`;

    await renameCategory(perTestCategoryName, newCategoryName);
    await assertCategoryDoesNotExist(perTestCategoryName);
    await assertCategoryExists(newCategoryName);
  });

  test("create a label", async () => {
    await setup(config);

    const labelName = `new-label-${config.tag}`;

    await assertLabelDoesNotExist(perTestCategoryName, labelName);

    await createLabel(perTestCategoryName, labelName);

    await assertLabelExists(perTestCategoryName, labelName);
  });

  test("delete a label", async () => {
    await setup(config);

    await deleteLabel(perTestCategoryName, perTestLabelName);
    await assertLabelDoesNotExist(perTestCategoryName, perTestLabelName);
  });

  test("rename a label", async () => {
    await setup(config);

    const newLabelName = "my-cool-new-label";

    await assertLabelDoesNotExist(perTestCategoryName, newLabelName);
    await renameLabel(perTestCategoryName, perTestLabelName, newLabelName);
    await assertLabelDoesNotExist(perTestCategoryName, perTestLabelName);
    await assertLabelExists(perTestCategoryName, newLabelName);
  });

  test("check cell count for a label loaded from file", async () => {
    await setup(config);

    const duplicateCategoryName = "duplicate";
    await duplicateCategory(duplicateCategoryName);

    await page.reload({ waitUntil: ["networkidle0", "domcontentloaded"] });

    const firstCategoryExpandIcon = await expect(page).toMatchElement(
      getTestClass("category-expand")
    );

    await firstCategoryExpandIcon.click();

    const expectedCategoryRow = await expect(page).toMatchElement(
      getTestClass("categorical-row")
    );
    const expectedLabelName = await getInnerText(
      expectedCategoryRow,
      "categorical-value"
    );
    const expectedLabelCount = await getInnerText(
      expectedCategoryRow,
      "categorical-value-count"
    );

    await expandCategory(duplicateCategoryName);

    const expectedCategory = await expect(page).toMatchElement(
      getTestClass("category")
    );

    const actualCategoryRow = await expect(expectedCategory).toMatchElement(
      getTestClass("categorical-row")
    );
    const actualLabelName = await getInnerText(
      actualCategoryRow,
      "categorical-value"
    );
    const actualLabelCount = await getInnerText(
      actualCategoryRow,
      "categorical-value-count"
    );

    expect(actualLabelName).toBe(expectedLabelName);
    expect(actualLabelCount).toBe(expectedLabelCount);

    async function getInnerText(element, className) {
      return element.$eval(getTestClass(className), (node) => node?.innerText);
    }
  });

  test("assign cells to a label", async () => {
    await setup(config);

    await expandCategory(perTestCategoryName);

    const lassoSelection = await calcDragCoordinates(
      "layout-graph",
      data.categoryLabel.lasso["coordinates-as-percent"]
    );

    await drag("layout-graph", lassoSelection.start, lassoSelection.end, true);
    await waitByID("lasso-element", { visible: true });
    await clickOn(`${perTestCategoryName}:${perTestLabelName}:see-actions`);
    await clickOn(
      `${perTestCategoryName}:${perTestLabelName}:add-current-selection-to-this-label`
    );

    const result = await waitByID(
      `categorical-value-count-${perTestCategoryName}-${perTestLabelName}`
    );

    expect(await result.evaluate((node) => node.innerText)).toBe(
      data.categoryLabel.newCount.bySubsetConfig[config.withSubset]
    );
  });

  test("undo/redo category creation", async () => {
    await setup(config);

    const categoryName = `category-created-undo-${config.tag}`;

    await assertCategoryDoesNotExist(categoryName);
    await createCategory(categoryName);
    await assertCategoryExists(categoryName);
    await clickOn("undo");
    await assertCategoryDoesNotExist(categoryName);
    await clickOn("redo");
    await assertCategoryExists(categoryName);
  });

  test("undo/redo category deletion", async () => {
    await setup(config);

    const categoryName = `category-deleted-undo-${config.tag}`;

    await createCategory(categoryName);
    await assertCategoryExists(categoryName);
    await deleteCategory(categoryName);
    await assertCategoryDoesNotExist(categoryName);
    await clickOn("undo");
    await assertCategoryExists(categoryName);
    await clickOn("redo");
    await assertCategoryDoesNotExist(categoryName);
  });

  test("undo/redo category rename", async () => {
    await setup(config);

    const newCategoryName = `category-renamed-undo-${config.tag}`;

    await assertCategoryDoesNotExist(newCategoryName);
    await renameCategory(perTestCategoryName, newCategoryName);
    await assertCategoryExists(newCategoryName);
    await assertCategoryDoesNotExist(perTestCategoryName);
    await clickOn("undo");
    await assertCategoryExists(perTestCategoryName);
    await assertCategoryDoesNotExist(newCategoryName);
    await clickOn("redo");
    await assertCategoryExists(newCategoryName);
    await assertCategoryDoesNotExist(perTestCategoryName);
  });

  test("undo/redo label creation", async () => {
    await setup(config);

    const labelName = `label-created-undo-${config.tag}`;

    await assertLabelDoesNotExist(perTestCategoryName, labelName);
    await createLabel(perTestCategoryName, labelName);
    await assertLabelExists(perTestCategoryName, labelName);
    await clickOn("undo");
    await assertLabelDoesNotExist(perTestCategoryName);
    await clickOn("redo");
    await assertLabelExists(perTestCategoryName, labelName);
  });

  test("undo/redo label deletion", async () => {
    await setup(config);

    await deleteLabel(perTestCategoryName, perTestLabelName);
    await assertLabelDoesNotExist(perTestCategoryName);
    await clickOn("undo");
    await assertLabelExists(perTestCategoryName, perTestLabelName);
    await clickOn("redo");
    await assertLabelDoesNotExist(perTestCategoryName);
  });

  test("undo/redo label rename", async () => {
    await setup(config);

    const newLabelName = `label-renamed-undo-${config.tag}`;

    await assertLabelDoesNotExist(perTestCategoryName, newLabelName);
    await renameLabel(perTestCategoryName, perTestLabelName, newLabelName);
    await assertLabelExists(perTestCategoryName, newLabelName);
    await assertLabelDoesNotExist(perTestCategoryName, perTestLabelName);
    await clickOn("undo");
    await assertLabelExists(perTestCategoryName, perTestLabelName);
    await assertLabelDoesNotExist(perTestCategoryName, newLabelName);
    await clickOn("redo");
    await assertLabelExists(perTestCategoryName, newLabelName);
    await assertLabelDoesNotExist(perTestCategoryName, perTestLabelName);
  });

  test("stacked bar graph renders", async () => {
    await setup(config);

    await expandCategory(perTestCategoryName);

    await clickOn(`colorby-louvain`);

    const labels = await getAllByClass("categorical-row");

    const result = await Promise.all(
      labels.map((label) => {
        return page.evaluate((element) => {
          return element.outerHTML;
        }, label);
      })
    );

    expect(result).toMatchSnapshot();
  });

  test("truncate midpoint whitespace", async () => {
    await setup(config);
    const newLabelName = "123 456";
    await renameLabel(perTestCategoryName, perTestLabelName, newLabelName);
    const value = await waitByID(
      `categorical-value-${perTestCategoryName}-${newLabelName}`
    );
    const result = await page.evaluate((elem) => elem.outerHTML, value);
    expect(result).toMatchSnapshot();
  });

  test("truncate single character", async () => {
    await setup(config);
    const newLabelName = "T";
    await renameLabel(perTestCategoryName, perTestLabelName, newLabelName);
    const value = await waitByID(
      `categorical-value-${perTestCategoryName}-${newLabelName}`
    );
    const result = await page.evaluate((elem) => elem.outerHTML, value);
    expect(result).toMatchSnapshot();
  });

  async function assertCategoryExists(categoryName) {
    const handle = await waitByID(`${categoryName}:category-label`);

    const result = await handle.evaluate((node) =>
      node.getAttribute("aria-label")
    );

    return expect(result).toBe(categoryName);
  }

  async function assertLabelExists(categoryName, labelName) {
    await expect(page).toMatchElement(
      getTestId(`${categoryName}:category-expand`)
    );

    await expandCategory(categoryName);

    const previous = await waitByID(
      `categorical-value-${categoryName}-${labelName}`
    );

    expect(
      await previous.evaluate((node) => node.getAttribute("aria-label"))
    ).toBe(labelName);
  }

  async function assertLabelDoesNotExist(categoryName, labelName) {
    await expandCategory(categoryName);
    const result = await page.$(
      `[data-testid='categorical-value-${categoryName}-${labelName}']`
    );
    expect(result).toBeNull();
  }
});
