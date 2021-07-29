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
  getOneElementInnerHTML,
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
  createGeneset,
  deleteGeneset,
  assertGenesetExists,
  assertGenesetDoesNotExist,
  getCellSetCount,
  expandGeneset,
  editGenesetName,
  addGeneToSet,
  assertGeneExistsInGeneset,
  removeGene,
  assertGeneDoesNotExist,
  expandGene,
  colorByGeneset,
  assertColorLegendLabel,
  colorByGene,
} from "./cellxgeneActions";

const data = datasets[DATASET];

const perTestCategoryName = "TEST-CATEGORY";
const perTestLabelName = "TEST-LABEL";

// geneset CRUD
const genesetToDeleteName = "geneset_to_delete";
const preExistingGenesetName = "fifth_dataset";
const meanExpressionBrushGenesetName = "second_gene_set";
const meanExpressionBrushCellsSelected = "557";
const subsetMeanExpressionBrushCellsSelected = "452";

// initial text, the text we type in, the result
const editableGenesetName = "geneset_to_edit";
const editText = "_111";
const newGenesetName = "geneset_to_edit_111";

// add gene to set
const geneToAddToSet = "RER1";
const setToAddGeneTo = "fill_this_geneset";

// remove gene from set
const geneToRemove = "SIK1";
const setToRemoveFrom = "empty_this_geneset";

// brush a gene
const geneToBrushAndColorBy = "SIK1";
const brushThisGeneGeneset = "brush_this_gene";
const geneBrushedCellCount = "109";
const subsetGeneBrushedCellCount = "96";

// eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
async function setup(config: any) {
  await goToPage(appUrlBase);

  if (config.categoricalAnno) {
    // setup the test fixtures
    await createCategory(perTestCategoryName);
    await createLabel(perTestCategoryName, perTestLabelName);
  }

  if (config.withSubset) {
    await subset({ x1: 0.1, y1: 0.1, x2: 0.8, y2: 0.8 });
  }

  await waitByClass("autosave-complete");
}

describe.each([
  { withSubset: true, tag: "subset" },
  { withSubset: false, tag: "whole" },
])("geneSET crud operations and interactions", (config) => {
  test("genesets load from csv", async () => {
    await setup(config);

    await assertGenesetExists(preExistingGenesetName);
  });
  test("brush on geneset mean", async () => {
    await setup(config);

    await expandGeneset(meanExpressionBrushGenesetName);

    const histBrushableAreaId = `histogram-${meanExpressionBrushGenesetName}-plot-brushable-area`;

    const coords = await calcDragCoordinates(histBrushableAreaId, {
      x1: 0.25,
      y1: 0.5,
      x2: 0.55,
      y2: 0.5,
    });

    await drag(histBrushableAreaId, coords.start, coords.end);

    const cellCount = await getCellSetCount(1);
    if (config.withSubset) {
      expect(cellCount).toBe(subsetMeanExpressionBrushCellsSelected);
    } else {
      expect(cellCount).toBe(meanExpressionBrushCellsSelected);
    }
  });
  test("color by mean expression", async () => {
    await setup(config);

    await colorByGeneset(meanExpressionBrushGenesetName);
    await assertColorLegendLabel(meanExpressionBrushGenesetName);
  });
  test("diffexp", async () => {
    if (config.withSubset) return;

    await setup(config);

    // set the two cell sets to b cells vs nk cells
    await expandCategory(`louvain`);
    await clickOn(`louvain:category-select`);
    await clickOn(`categorical-value-select-louvain-B cells`);
    await clickOn(`cellset-button-1`);
    await clickOn(`categorical-value-select-louvain-B cells`);
    await clickOn(`categorical-value-select-louvain-NK cells`);
    await clickOn(`cellset-button-2`);

    // run diffexp
    await clickOn(`diffexp-button`);
    await waitByClass("pop-1-geneset-expand");
    await expect(page).toClick(getTestClass("pop-1-geneset-expand"));

    await page.waitForFunction(
      // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
      (selector: any) => !document.querySelector(selector),
      {},
      getTestClass("gene-loading-spinner")
    );

    let genesHTML = await getOneElementInnerHTML(
      getTestClass("gene-set-genes")
    );

    expect(genesHTML).toMatchSnapshot();

    await expect(page).toClick(getTestClass("pop-1-geneset-expand"));
    await expect(page).toClick(getTestClass("pop-2-geneset-expand"));

    await page.waitForFunction(
      // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
      (selector: any) => !document.querySelector(selector),
      {},
      getTestClass("gene-loading-spinner")
    );

    genesHTML = await getOneElementInnerHTML(getTestClass("gene-set-genes"));

    expect(genesHTML).toMatchSnapshot();
  });
  test("create a new geneset", async () => {
    if (config.withSubset) return;

    await setup(config);

    const genesetName = `test-geneset-foo-123`;
    await assertGenesetDoesNotExist(genesetName);
    await createGeneset(genesetName);
    /* note: as of June 2021, the aria label is in the truncate component which clones the element */
    await assertGenesetExists(genesetName);
  });
  test("edit geneset name", async () => {
    await setup(config);

    await editGenesetName(editableGenesetName, editText);
    await assertGenesetExists(newGenesetName);
  });
  test("delete a geneset", async () => {
    if (config.withSubset) return;

    await setup(config);

    await deleteGeneset(genesetToDeleteName);
  });
});

describe.each([
  { withSubset: true, tag: "subset" },
  { withSubset: false, tag: "whole" },
])("GENE crud operations and interactions", (config) => {
  test("add a gene to geneset", async () => {
    await setup(config);

    await addGeneToSet(setToAddGeneTo, geneToAddToSet);
    await expandGeneset(setToAddGeneTo);
    await assertGeneExistsInGeneset(geneToAddToSet);
  });
  test("expand gene and brush", async () => {
    await setup(config);

    await expandGeneset(brushThisGeneGeneset);
    await expandGene(geneToBrushAndColorBy);
    const histBrushableAreaId = `histogram-${geneToBrushAndColorBy}-plot-brushable-area`;

    const coords = await calcDragCoordinates(histBrushableAreaId, {
      x1: 0.25,
      y1: 0.5,
      x2: 0.55,
      y2: 0.5,
    });
    await drag(histBrushableAreaId, coords.start, coords.end);
    const cellCount = await getCellSetCount(1);
    if (config.withSubset) {
      expect(cellCount).toBe(subsetGeneBrushedCellCount);
    } else {
      expect(cellCount).toBe(geneBrushedCellCount);
    }
  });
  test("color by gene in geneset", async () => {
    await setup(config);

    await expandGeneset(meanExpressionBrushGenesetName);

    await colorByGene(geneToBrushAndColorBy);
    await assertColorLegendLabel(geneToBrushAndColorBy);
  });
  test("delete gene from geneset", async () => {
    // We've already deleted the gene
    if (config.withSubset) return;

    await setup(config);

    await expandGeneset(setToRemoveFrom);
    await removeGene(geneToRemove);
    await assertGeneDoesNotExist(geneToRemove);
  });
});

describe.each([
  { withSubset: true, tag: "subset", categoricalAnno: true },
  { withSubset: false, tag: "whole", categoricalAnno: true },
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

    // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
    async function getInnerText(element: any, className: any) {
      return element.$eval(
        getTestClass(className),
        // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
        (node: any) => node?.innerText
      );
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

    // @ts-expect-error ts-migrate(2531) FIXME: Object is possibly 'null'.
    expect(await result.evaluate((node) => node.innerText)).toBe(
      // @ts-expect-error ts-migrate(2538) FIXME: Type 'boolean' cannot be used as an index type.
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
    // @ts-expect-error ts-migrate(2554) FIXME: Expected 2 arguments, but got 1.
    await assertLabelDoesNotExist(perTestCategoryName);
    await clickOn("redo");
    await assertLabelExists(perTestCategoryName, labelName);
  });

  test("undo/redo label deletion", async () => {
    await setup(config);

    await deleteLabel(perTestCategoryName, perTestLabelName);
    // @ts-expect-error ts-migrate(2554) FIXME: Expected 2 arguments, but got 1.
    await assertLabelDoesNotExist(perTestCategoryName);
    await clickOn("undo");
    await assertLabelExists(perTestCategoryName, perTestLabelName);
    await clickOn("redo");
    // @ts-expect-error ts-migrate(2554) FIXME: Expected 2 arguments, but got 1.
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

  // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
  async function assertCategoryExists(categoryName: any) {
    const handle = await waitByID(`${categoryName}:category-label`);

    // @ts-expect-error ts-migrate(2531) FIXME: Object is possibly 'null'.
    const result = await handle.evaluate((node) =>
      node.getAttribute("aria-label")
    );

    return expect(result).toBe(categoryName);
  }

  // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
  async function assertLabelExists(categoryName: any, labelName: any) {
    await expect(page).toMatchElement(
      getTestId(`${categoryName}:category-expand`)
    );

    await expandCategory(categoryName);

    const previous = await waitByID(
      `categorical-value-${categoryName}-${labelName}`
    );

    expect(
      // @ts-expect-error ts-migrate(2531) FIXME: Object is possibly 'null'.
      await previous.evaluate((node) => node.getAttribute("aria-label"))
    ).toBe(labelName);
  }

  // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
  async function assertLabelDoesNotExist(categoryName: any, labelName: any) {
    await expandCategory(categoryName);
    const result = await page.$(
      `[data-testid='categorical-value-${categoryName}-${labelName}']`
    );
    expect(result).toBeNull();
  }
});
