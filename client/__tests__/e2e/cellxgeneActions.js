/* eslint-disable no-await-in-loop -- await in loop is needed to emulate sequential user actions  */
import { strict as assert } from "assert";
import {
  clearInputAndTypeInto,
  clickOn,
  getAllByClass,
  getOneElementInnerText,
  typeInto,
  waitByID,
  waitByClass,
  waitForAllByIds,
  clickOnUntil,
  getTestClass,
  getTestId,
  isElementPresent,
  goToPage,
} from "./puppeteerUtils";

import { appUrlBase, TEST_EMAIL, TEST_PASSWORD } from "./config";

export async function drag(testId, start, end, lasso = false) {
  const layout = await waitByID(testId);
  const elBox = await layout.boxModel();
  const x1 = elBox.content[0].x + start.x;
  const x2 = elBox.content[0].x + end.x;
  const y1 = elBox.content[0].y + start.y;
  const y2 = elBox.content[0].y + end.y;
  await page.mouse.move(x1, y1);
  await page.mouse.down();
  if (lasso) {
    await page.mouse.move(x2, y1);
    await page.mouse.move(x2, y2);
    await page.mouse.move(x1, y2);
    await page.mouse.move(x1, y1);
  } else {
    await page.mouse.move(x2, y2);
  }
  await page.mouse.up();
}

export async function clickOnCoordinate(testId, coord) {
  const layout = await expect(page).toMatchElement(getTestId(testId));
  const elBox = await layout.boxModel();

  if (!elBox) {
    throw Error("Layout's boxModel is not available!");
  }

  const x = elBox.content[0].x + coord.x;
  const y = elBox.content[0].y + coord.y;
  await page.mouse.click(x, y);
}

export async function getAllHistograms(testclass, testIds) {
  const histTestIds = testIds.map((tid) => `histogram-${tid}`);

  // these load asynchronously, so we need to wait for each histogram individually,
  // and they may be quite slow in some cases.
  await waitForAllByIds(histTestIds, { timeout: 4 * 60 * 1000 });

  const allHistograms = await getAllByClass(testclass);

  const testIDs = await Promise.all(
    allHistograms.map((hist) => page.evaluate((elem) => elem.dataset.testid, hist))
  );

  return testIDs.map((id) => id.replace(/^histogram-/, ""));
}

export async function getAllCategoriesAndCounts(category) {
  // these load asynchronously, so we have to wait for the specific category.
  await waitByID(`category-${category}`);

  return page.$$eval(
    `[data-testid="category-${category}"] [data-testclass='categorical-row']`,
    (rows) =>
      Object.fromEntries(
        rows.map((row) => {
          const cat = row
            .querySelector("[data-testclass='categorical-value']")
            .getAttribute("aria-label");

          const count = row.querySelector(
            "[data-testclass='categorical-value-count']"
          ).innerText;

          return [cat, count];
        })
      )
  );
}

export async function getCellSetCount(num) {
  await clickOn(`cellset-button-${num}`);
  return getOneElementInnerText(`[data-testid='cellset-count-${num}']`);
}

export async function resetCategory(category) {
  const checkboxId = `${category}:category-select`;
  await waitByID(checkboxId);
  const checkedPseudoclass = await page.$eval(
    `[data-testid='${checkboxId}']`,
    (el) => el.matches(":checked")
  );
  if (!checkedPseudoclass) await clickOn(checkboxId);

  const categoryRow = await waitByID(`${category}:category-expand`);

  const isExpanded = await categoryRow.$(
    "[data-testclass='category-expand-is-expanded']"
  );

  if (isExpanded) await clickOn(`${category}:category-expand`);
}

export async function calcCoordinate(testId, xAsPercent, yAsPercent) {
  const el = await waitByID(testId);
  const size = await el.boxModel();
  return {
    x: Math.floor(size.width * xAsPercent),
    y: Math.floor(size.height * yAsPercent),
  };
}

export async function calcDragCoordinates(testId, coordinateAsPercent) {
  return {
    start: await calcCoordinate(
      testId,
      coordinateAsPercent.x1,
      coordinateAsPercent.y1
    ),
    end: await calcCoordinate(
      testId,
      coordinateAsPercent.x2,
      coordinateAsPercent.y2
    ),
  };
}

export async function selectCategory(category, values, reset = true) {
  if (reset) await resetCategory(category);

  await clickOn(`${category}:category-expand`);
  await clickOn(`${category}:category-select`);

  for (const value of values) {
    await clickOn(`categorical-value-select-${category}-${value}`);
  }
}

export async function expandCategory(category) {
  const expand = await waitByID(`${category}:category-expand`);
  const notExpanded = await expand.$(
    "[data-testclass='category-expand-is-not-expanded']"
  );
  if (notExpanded) await clickOn(`${category}:category-expand`);
}

export async function clip(min = 0, max = 100) {
  await clickOn("visualization-settings");
  await clearInputAndTypeInto("clip-min-input", min);
  await clearInputAndTypeInto("clip-max-input", max);
  await clickOn("clip-commit");
}

export async function createCategory(categoryName) {
  await clickOnUntil("open-annotation-dialog", async () => {
    await expect(page).toMatchElement(getTestId("new-category-name"));
  });

  await typeInto("new-category-name", categoryName);
  await clickOn("submit-category");
}

/* 

  GENESET 

*/

export async function colorByGeneset(genesetName) {
  await clickOn(`${genesetName}:colorby-entire-geneset`);
}

export async function colorByGene(gene) {
  await clickOn(`colorby-${gene}`);
}

export async function assertColorLegendLabel(label) {
  const handle = await waitByID("continuous_legend_color_by_label");

  const result = await handle.evaluate((node) => node.getAttribute("aria-label"));

  return expect(result).toBe(label);
}

export async function expandGeneset(genesetName) {
  const expand = await waitByID(`${genesetName}:geneset-expand`);
  const notExpanded = await expand.$(
    "[data-testclass='geneset-expand-is-not-expanded']"
  );
  if (notExpanded) await clickOn(`${genesetName}:geneset-expand`);
}

export async function createGeneset(genesetName) {
  await clickOnUntil("open-create-geneset-dialog", async () => {
    await expect(page).toMatchElement(getTestId("create-geneset-input"));
  });

  await typeInto("create-geneset-input", genesetName);
  await clickOn("submit-geneset");
  await waitByClass("autosave-complete");
}

export async function editGenesetName(genesetName, editText) {
  const editButton = `${genesetName}:edit-genesetName-mode`;
  const submitButton = `${genesetName}:submit-geneset`;
  await clickOnUntil(`${genesetName}:see-actions`, async () => {
    await expect(page).toMatchElement(getTestId(editButton));
  });
  await clickOn(editButton);
  await typeInto("rename-geneset-modal", editText);
  await clickOn(submitButton);
}

export async function deleteGeneset(genesetName) {
  const targetId = `${genesetName}:delete-geneset`;

  await clickOnUntil(`${genesetName}:see-actions`, async () => {
    await expect(page).toMatchElement(getTestId(targetId));
  });

  await clickOn(targetId);

  await assertGenesetDoesNotExist(genesetName);
  await waitByClass("autosave-complete");
}

export async function assertGenesetDoesNotExist(genesetName) {
  const result = await isElementPresent(
    getTestId(`${genesetName}:geneset-name`)
  );
  await expect(result).toBe(false);
}

export async function assertGenesetExists(genesetName) {
  const handle = await waitByID(`${genesetName}:geneset-name`);

  const result = await handle.evaluate((node) => node.getAttribute("aria-label"));

  return expect(result).toBe(genesetName);
}

/* 

  GENE

*/

export async function addGeneToSet(genesetName, geneToAddToSet) {
  const submitButton = `${genesetName}:submit-gene`;

  await clickOn(`${genesetName}:add-new-gene-to-geneset`);
  await typeInto("add-genes", geneToAddToSet);
  await clickOn(submitButton);
}

export async function removeGene(geneSymbol) {
  const targetId = `delete-from-geneset:${geneSymbol}`;

  await clickOn(targetId);

  await waitByClass("autosave-complete");
}

export async function assertGeneExistsInGeneset(geneSymbol) {
  const handle = await waitByID(`${geneSymbol}:gene-label`);

  const result = await handle.evaluate((node) => node.getAttribute("aria-label"));

  return expect(result).toBe(geneSymbol);
}

export async function assertGeneDoesNotExist(geneSymbol) {
  const result = await isElementPresent(getTestId(`${geneSymbol}:gene-label`));

  await expect(result).toBe(false);
}

export async function expandGene(geneSymbol) {
  await clickOn(`maximize-${geneSymbol}`);
}

/* 

  CATEGORY 

*/

export async function duplicateCategory(categoryName) {
  await clickOn("open-annotation-dialog");

  await typeInto("new-category-name", categoryName);

  const dropdownOptionClass = "duplicate-category-dropdown-option";

  await clickOnUntil("duplicate-category-dropdown", async () => {
    await expect(page).toMatchElement(getTestClass(dropdownOptionClass));
  });

  const option = await expect(page).toMatchElement(
    getTestClass(dropdownOptionClass)
  );

  await option.click();

  await clickOnUntil("submit-category", async () => {
    await expect(page).toMatchElement(
      getTestId(`${categoryName}:category-expand`)
    );
  });

  await waitByClass("autosave-complete");
}

export async function renameCategory(oldCategoryName, newCategoryName) {
  await clickOn(`${oldCategoryName}:see-actions`);
  await clickOn(`${oldCategoryName}:edit-category-mode`);
  await clearInputAndTypeInto(
    `${oldCategoryName}:edit-category-name-text`,
    newCategoryName
  );
  await clickOn(`${oldCategoryName}:submit-category-edit`);
}

export async function deleteCategory(categoryName) {
  const targetId = `${categoryName}:delete-category`;

  await clickOnUntil(`${categoryName}:see-actions`, async () => {
    await expect(page).toMatchElement(getTestId(targetId));
  });

  await clickOn(targetId);

  await assertCategoryDoesNotExist();
}

export async function createLabel(categoryName, labelName) {
  /**
   * (thuang): This explicit wait is needed, since currently showing
   * the modal again quickly after the previous action dismissing the
   * modal will persist the input value from the previous action.
   *
   * To reproduce:
   * 1. Click on the plus sign to show the modal to add a new label to the category
   * 2. Type `123` in the input box
   * 3. Hover over your mouse over the plus sign and double click to quickly dismiss and
   * invoke the modal again
   * 4. You will see `123` is persisted in the input box
   * 5. Expected behavior is to get an empty input box
   */
  await page.waitForTimeout(500);

  await clickOn(`${categoryName}:see-actions`);

  await clickOn(`${categoryName}:add-new-label-to-category`);

  await typeInto(`${categoryName}:new-label-name`, labelName);

  await clickOn(`${categoryName}:submit-label`);
}

export async function deleteLabel(categoryName, labelName) {
  await expandCategory(categoryName);
  await clickOn(`${categoryName}:${labelName}:see-actions`);
  await clickOn(`${categoryName}:${labelName}:delete-label`);
}

export async function renameLabel(categoryName, oldLabelName, newLabelName) {
  await expandCategory(categoryName);
  await clickOn(`${categoryName}:${oldLabelName}:see-actions`);
  await clickOn(`${categoryName}:${oldLabelName}:edit-label`);
  await clearInputAndTypeInto(
    `${categoryName}:${oldLabelName}:edit-label-name`,
    newLabelName
  );
  await clickOn(`${categoryName}:${oldLabelName}:submit-label-edit`);
}

export async function addGeneToSearch(geneName) {
  await typeInto("gene-search", geneName);
  await page.keyboard.press("Enter");
  await page.waitForSelector(`[data-testid='histogram-${geneName}']`);
}

export async function subset(coordinatesAsPercent) {
  // In order to deselect the selection after the subset, make sure we have some clear part
  // of the scatterplot we can click on
  assert(coordinatesAsPercent.x2 < 0.99 || coordinatesAsPercent.y2 < 0.99);
  const lassoSelection = await calcDragCoordinates(
    "layout-graph",
    coordinatesAsPercent
  );
  await drag("layout-graph", lassoSelection.start, lassoSelection.end, true);
  await clickOn("subset-button");
  const clearCoordinate = await calcCoordinate("layout-graph", 0.5, 0.99);
  await clickOnCoordinate("layout-graph", clearCoordinate);
}

export async function setSellSet(cellSet, cellSetNum) {
  const selections = cellSet.filter((sel) => sel.kind === "categorical");

  for (const selection of selections) {
    await selectCategory(selection.metadata, selection.values, true);
  }

  await getCellSetCount(cellSetNum);
}

export async function runDiffExp(cellSet1, cellSet2) {
  await setSellSet(cellSet1, 1);
  await setSellSet(cellSet2, 2);
  await clickOn("diffexp-button");
}

export async function bulkAddGenes(geneNames) {
  await clickOn("section-bulk-add");
  await typeInto("input-bulk-add", geneNames.join(","));
  await page.keyboard.press("Enter");
}

export async function assertCategoryDoesNotExist(categoryName) {
  const result = await isElementPresent(
    getTestId(`${categoryName}:category-label`)
  );

  await expect(result).toBe(false);
}

export async function login() {
  await goToPage(appUrlBase);

  await clickOn("log-in");

  // (thuang): Auth0 form is unstable and unsafe for input until verified
  await waitUntilFormFieldStable('[name="email"]');

  await expect(page).toFillForm("form", {
    email: TEST_EMAIL,
    password: TEST_PASSWORD,
  });

  await Promise.all([
    page.waitForNavigation({ waitUntil: "networkidle0" }),
    expect(page).toClick('[name="submit"]'),
  ]);

  expect(page.url()).toContain(appUrlBase);
}

export async function logout() {
  await clickOnUntil("user-info", async () => {
    await waitByID("log-out");
    await Promise.all([
      page.waitForNavigation({ waitUntil: "networkidle0" }),
      clickOn("log-out"),
    ]);
  });

  await waitByID("log-in");
}

async function waitUntilFormFieldStable(selector) {
  const MAX_RETRY = 10;
  const WAIT_FOR_MS = 200;

  const EXPECTED_VALUE = "aaa";

  let retry = 0;

  while (retry < MAX_RETRY) {
    try {
      await expect(page).toFill(selector, EXPECTED_VALUE);

      const fieldHandle = await expect(page).toMatchElement(selector);

      const fieldValue = await page.evaluate(
        (input) => input.value,
        fieldHandle
      );

      expect(fieldValue).toBe(EXPECTED_VALUE);

      break;
    } catch (error) {
      retry += 1;

      await page.waitForTimeout(WAIT_FOR_MS);
    }
  }

  if (retry === MAX_RETRY) {
    throw Error("clickOnUntil() assertion failed!");
  }
}
/* eslint-enable no-await-in-loop -- await in loop is needed to emulate sequential user actions */
