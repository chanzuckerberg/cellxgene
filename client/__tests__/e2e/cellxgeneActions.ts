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

// eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
export async function drag(testId: any, start: any, end: any, lasso = false) {
  const layout = await waitByID(testId);
  // @ts-expect-error ts-migrate(2531) FIXME: Object is possibly 'null'.
  const elBox = await layout.boxModel();
  // @ts-expect-error ts-migrate(2531) FIXME: Object is possibly 'null'.
  const x1 = elBox.content[0].x + start.x;
  // @ts-expect-error ts-migrate(2531) FIXME: Object is possibly 'null'.
  const x2 = elBox.content[0].x + end.x;
  // @ts-expect-error ts-migrate(2531) FIXME: Object is possibly 'null'.
  const y1 = elBox.content[0].y + start.y;
  // @ts-expect-error ts-migrate(2531) FIXME: Object is possibly 'null'.
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

// eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
export async function clickOnCoordinate(testId: any, coord: any) {
  const layout = await expect(page).toMatchElement(getTestId(testId));
  const elBox = await layout.boxModel();

  if (!elBox) {
    throw Error("Layout's boxModel is not available!");
  }

  const x = elBox.content[0].x + coord.x;
  const y = elBox.content[0].y + coord.y;
  await page.mouse.click(x, y);
}

// eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
export async function getAllHistograms(testclass: any, testIds: any) {
  // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
  const histTestIds = testIds.map((tid: any) => `histogram-${tid}`);

  // these load asynchronously, so we need to wait for each histogram individually,
  // and they may be quite slow in some cases.
  // @ts-expect-error ts-migrate(2554) FIXME: Expected 1 arguments, but got 2.
  await waitForAllByIds(histTestIds, { timeout: 4 * 60 * 1000 });

  const allHistograms = await getAllByClass(testclass);

  const testIDs = await Promise.all(
    allHistograms.map((hist) => {
      return page.evaluate((elem) => {
        return elem.dataset.testid;
      }, hist);
    })
  );

  return testIDs.map((id) => id.replace(/^histogram-/, ""));
}

// eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
export async function getAllCategoriesAndCounts(category: any) {
  // these load asynchronously, so we have to wait for the specific category.
  await waitByID(`category-${category}`);

  return page.$$eval(
    `[data-testid="category-${category}"] [data-testclass='categorical-row']`,
    (rows) =>
      Object.fromEntries(
        rows.map((row) => {
          // @ts-expect-error ts-migrate(2531) FIXME: Object is possibly 'null'.
          const cat = row
            .querySelector("[data-testclass='categorical-value']")
            .getAttribute("aria-label");

          const count = (row.querySelector(
            "[data-testclass='categorical-value-count']"
            // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
          ) as any).innerText;

          return [cat, count];
        })
      )
  );
}

// eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
export async function getCellSetCount(num: any) {
  await clickOn(`cellset-button-${num}`);
  return getOneElementInnerText(`[data-testid='cellset-count-${num}']`);
}

// eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
export async function resetCategory(category: any) {
  const checkboxId = `${category}:category-select`;
  await waitByID(checkboxId);
  const checkedPseudoclass = await page.$eval(
    `[data-testid='${checkboxId}']`,
    (el) => el.matches(":checked")
  );
  if (!checkedPseudoclass) await clickOn(checkboxId);

  const categoryRow = await waitByID(`${category}:category-expand`);

  // @ts-expect-error ts-migrate(2531) FIXME: Object is possibly 'null'.
  const isExpanded = await categoryRow.$(
    "[data-testclass='category-expand-is-expanded']"
  );

  if (isExpanded) await clickOn(`${category}:category-expand`);
}

// eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types --- FIXME: disabled temporarily on migrate to TS.
export // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
async function calcCoordinate(testId: any, xAsPercent: any, yAsPercent: any) {
  const el = await waitByID(testId);
  // @ts-expect-error ts-migrate(2531) FIXME: Object is possibly 'null'.
  const size = await el.boxModel();
  return {
    // @ts-expect-error ts-migrate(2531) FIXME: Object is possibly 'null'.
    x: Math.floor(size.width * xAsPercent),
    // @ts-expect-error ts-migrate(2531) FIXME: Object is possibly 'null'.
    y: Math.floor(size.height * yAsPercent),
  };
}

// eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types --- FIXME: disabled temporarily on migrate to TS.
export // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
async function calcDragCoordinates(testId: any, coordinateAsPercent: any) {
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

// eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
export async function selectCategory(category: any, values: any, reset = true) {
  if (reset) await resetCategory(category);

  await clickOn(`${category}:category-expand`);
  await clickOn(`${category}:category-select`);

  for (const value of values) {
    await clickOn(`categorical-value-select-${category}-${value}`);
  }
}

// eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
export async function expandCategory(category: any) {
  const expand = await waitByID(`${category}:category-expand`);
  // @ts-expect-error ts-migrate(2531) FIXME: Object is possibly 'null'.
  const notExpanded = await expand.$(
    "[data-testclass='category-expand-is-not-expanded']"
  );
  if (notExpanded) await clickOn(`${category}:category-expand`);
}

// eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types --- FIXME: disabled temporarily on migrate to TS.
export async function clip(min = 0, max = 100) {
  await clickOn("visualization-settings");
  await clearInputAndTypeInto("clip-min-input", min);
  await clearInputAndTypeInto("clip-max-input", max);
  await clickOn("clip-commit");
}

// eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
export async function createCategory(categoryName: any) {
  await clickOnUntil("open-annotation-dialog", async () => {
    await expect(page).toMatchElement(getTestId("new-category-name"));
  });

  await typeInto("new-category-name", categoryName);
  await clickOn("submit-category");
}

/* 

  GENESET 

*/

// eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
export async function colorByGeneset(genesetName: any) {
  await clickOn(`${genesetName}:colorby-entire-geneset`);
}

// eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
export async function colorByGene(gene: any) {
  await clickOn(`colorby-${gene}`);
}

// eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
export async function assertColorLegendLabel(label: any) {
  const handle = await waitByID("continuous_legend_color_by_label");

  // @ts-expect-error ts-migrate(2531) FIXME: Object is possibly 'null'.
  const result = await handle.evaluate((node) => {
    return node.getAttribute("aria-label");
  });

  return expect(result).toBe(label);
}

// eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
export async function expandGeneset(genesetName: any) {
  const expand = await waitByID(`${genesetName}:geneset-expand`);
  // @ts-expect-error ts-migrate(2531) FIXME: Object is possibly 'null'.
  const notExpanded = await expand.$(
    "[data-testclass='geneset-expand-is-not-expanded']"
  );
  if (notExpanded) await clickOn(`${genesetName}:geneset-expand`);
}

// eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
export async function createGeneset(genesetName: any) {
  await clickOnUntil("open-create-geneset-dialog", async () => {
    await expect(page).toMatchElement(getTestId("create-geneset-input"));
  });

  await typeInto("create-geneset-input", genesetName);
  await clickOn("submit-geneset");
  await waitByClass("autosave-complete");
}

// eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
export async function editGenesetName(genesetName: any, editText: any) {
  const editButton = `${genesetName}:edit-genesetName-mode`;
  const submitButton = `${genesetName}:submit-geneset`;
  await clickOnUntil(`${genesetName}:see-actions`, async () => {
    await expect(page).toMatchElement(getTestId(editButton));
  });
  await clickOn(editButton);
  await typeInto("rename-geneset-modal", editText);
  await clickOn(submitButton);
}

// eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
export async function deleteGeneset(genesetName: any) {
  const targetId = `${genesetName}:delete-geneset`;

  await clickOnUntil(`${genesetName}:see-actions`, async () => {
    await expect(page).toMatchElement(getTestId(targetId));
  });

  await clickOn(targetId);

  await assertGenesetDoesNotExist(genesetName);
  await waitByClass("autosave-complete");
}

// eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
export async function assertGenesetDoesNotExist(genesetName: any) {
  // @ts-expect-error ts-migrate(2554) FIXME: Expected 2 arguments, but got 1.
  const result = await isElementPresent(
    getTestId(`${genesetName}:geneset-name`)
  );
  await expect(result).toBe(false);
}

// eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
export async function assertGenesetExists(genesetName: any) {
  const handle = await waitByID(`${genesetName}:geneset-name`);

  // @ts-expect-error ts-migrate(2531) FIXME: Object is possibly 'null'.
  const result = await handle.evaluate((node) => {
    return node.getAttribute("aria-label");
  });

  return expect(result).toBe(genesetName);
}

/* 

  GENE

*/

// eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
export async function addGeneToSet(genesetName: any, geneToAddToSet: any) {
  const submitButton = `${genesetName}:submit-gene`;

  await clickOn(`${genesetName}:add-new-gene-to-geneset`);
  await typeInto("add-genes", geneToAddToSet);
  await clickOn(submitButton);
}

// eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
export async function removeGene(geneSymbol: any) {
  const targetId = `delete-from-geneset:${geneSymbol}`;

  await clickOn(targetId);

  await waitByClass("autosave-complete");
}

// eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
export async function assertGeneExistsInGeneset(geneSymbol: any) {
  const handle = await waitByID(`${geneSymbol}:gene-label`);

  // @ts-expect-error ts-migrate(2531) FIXME: Object is possibly 'null'.
  const result = await handle.evaluate((node) => {
    return node.getAttribute("aria-label");
  });

  return expect(result).toBe(geneSymbol);
}

// eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
export async function assertGeneDoesNotExist(geneSymbol: any) {
  // @ts-expect-error ts-migrate(2554) FIXME: Expected 2 arguments, but got 1.
  const result = await isElementPresent(getTestId(`${geneSymbol}:gene-label`));

  await expect(result).toBe(false);
}

// eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
export async function expandGene(geneSymbol: any) {
  await clickOn(`maximize-${geneSymbol}`);
}

/* 

  CATEGORY 

*/

// eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
export async function duplicateCategory(categoryName: any) {
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

// eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types --- FIXME: disabled temporarily on migrate to TS.
export // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
async function renameCategory(oldCategoryName: any, newCategoryName: any) {
  await clickOn(`${oldCategoryName}:see-actions`);
  await clickOn(`${oldCategoryName}:edit-category-mode`);
  await clearInputAndTypeInto(
    `${oldCategoryName}:edit-category-name-text`,
    newCategoryName
  );
  await clickOn(`${oldCategoryName}:submit-category-edit`);
}

// eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
export async function deleteCategory(categoryName: any) {
  const targetId = `${categoryName}:delete-category`;

  await clickOnUntil(`${categoryName}:see-actions`, async () => {
    await expect(page).toMatchElement(getTestId(targetId));
  });

  await clickOn(targetId);

  // @ts-expect-error ts-migrate(2554) FIXME: Expected 1 arguments, but got 0.
  await assertCategoryDoesNotExist();
}

// eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
export async function createLabel(categoryName: any, labelName: any) {
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

// eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
export async function deleteLabel(categoryName: any, labelName: any) {
  await expandCategory(categoryName);
  await clickOn(`${categoryName}:${labelName}:see-actions`);
  await clickOn(`${categoryName}:${labelName}:delete-label`);
}

// eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types --- FIXME: disabled temporarily on migrate to TS.
export async function renameLabel(
  // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
  categoryName: any,
  // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
  oldLabelName: any,
  // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
  newLabelName: any
) {
  await expandCategory(categoryName);
  await clickOn(`${categoryName}:${oldLabelName}:see-actions`);
  await clickOn(`${categoryName}:${oldLabelName}:edit-label`);
  await clearInputAndTypeInto(
    `${categoryName}:${oldLabelName}:edit-label-name`,
    newLabelName
  );
  await clickOn(`${categoryName}:${oldLabelName}:submit-label-edit`);
}

// eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
export async function addGeneToSearch(geneName: any) {
  await typeInto("gene-search", geneName);
  await page.keyboard.press("Enter");
  await page.waitForSelector(`[data-testid='histogram-${geneName}']`);
}

// eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
export async function subset(coordinatesAsPercent: any) {
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

// eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
export async function setSellSet(cellSet: any, cellSetNum: any) {
  // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
  const selections = cellSet.filter((sel: any) => sel.kind === "categorical");

  for (const selection of selections) {
    await selectCategory(selection.metadata, selection.values, true);
  }

  await getCellSetCount(cellSetNum);
}

// eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
export async function runDiffExp(cellSet1: any, cellSet2: any) {
  await setSellSet(cellSet1, 1);
  await setSellSet(cellSet2, 2);
  await clickOn("diffexp-button");
}

// eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
export async function bulkAddGenes(geneNames: any) {
  await clickOn("section-bulk-add");
  await typeInto("input-bulk-add", geneNames.join(","));
  await page.keyboard.press("Enter");
}

// eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
export async function assertCategoryDoesNotExist(categoryName: any) {
  // @ts-expect-error ts-migrate(2554) FIXME: Expected 2 arguments, but got 1.
  const result = await isElementPresent(
    getTestId(`${categoryName}:category-label`)
  );

  await expect(result).toBe(false);
}

// eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types --- FIXME: disabled temporarily on migrate to TS.
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

// eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types --- FIXME: disabled temporarily on migrate to TS.
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

// eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
async function waitUntilFormFieldStable(selector: any) {
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
