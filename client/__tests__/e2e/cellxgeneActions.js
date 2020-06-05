import { strict as assert } from "assert";

export const cellxgeneActions = (page, utils) => ({
  async drag(testId, start, end, lasso = false) {
    const layout = await utils.waitByID(testId);
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
  },

  async clickOnCoordinate(testId, coord) {
    const layout = await utils.waitByID(testId);
    const elBox = await layout.boxModel();
    const x = elBox.content[0].x + coord.x;
    const y = elBox.content[0].y + coord.y;
    await page.mouse.click(x, y);
  },

  async getAllHistograms(testclass, testIds) {
    const histTestIds = testIds.map((tid) => `histogram-${tid}`);
    // these load asynchronously, so we need to wait for each histogram individually
    await utils.waitForAllByIds(histTestIds);
    const allHistograms = await utils.getAllByClass(testclass);
    return allHistograms.map((hist) => hist.replace(/^histogram-/, ""));
  },

  async getAllCategoriesAndCounts(category) {
    await utils.waitByClass("categorical-row");
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
  },

  async cellSet(num) {
    await utils.clickOn(`cellset-button-${num}`);
    return utils.getOneElementInnerText(`[data-testid='cellset-count-${num}']`);
  },

  async resetCategory(category) {
    const checkboxId = `${category}:category-select`;
    await utils.waitByID(checkboxId);
    const checkedPseudoclass = await page.$eval(
      `[data-testid='${checkboxId}']`,
      (el) => el.matches(":checked")
    );
    if (!checkedPseudoclass) await utils.clickOn(checkboxId);
    try {
      const categoryRow = await utils.waitByID(`${category}:category-expand`);
      const isExpanded = await categoryRow.$(
        "[data-testclass='category-expand-is-expanded']"
      );
      if (isExpanded) await utils.clickOn(`${category}:category-expand`);
    } catch {}
  },

  async calcCoordinate(testId, xAsPercent, yAsPercent) {
    const el = await utils.waitByID(testId);
    const size = await el.boxModel();
    return {
      x: Math.floor(size.width * xAsPercent),
      y: Math.floor(size.height * yAsPercent),
    };
  },

  async calcDragCoordinates(testId, coordinateAsPercent) {
    return {
      start: await this.calcCoordinate(
        testId,
        coordinateAsPercent.x1,
        coordinateAsPercent.y1
      ),
      end: await this.calcCoordinate(
        testId,
        coordinateAsPercent.x2,
        coordinateAsPercent.y2
      ),
    };
  },

  async selectCategory(category, values, reset = true) {
    if (reset) await this.resetCategory(category);
    await utils.clickOn(`${category}:category-expand`);
    await utils.clickOn(`${category}:category-select`);
    for (const val of values) {
      await utils.clickOn(`categorical-value-select-${category}-${val}`);
    }
  },

  async expandCategory(category) {
    const expand = await utils.waitByID(`${category}:category-expand`);
    const notExpanded = await expand.$(
      "[data-testclass='category-expand-is-not-expanded']"
    );
    if (notExpanded) await utils.clickOn(`${category}:category-expand`);
  },

  async clip(min = 0, max = 100) {
    await utils.clickOn("visualization-settings");
    await utils.clearInputAndTypeInto("clip-min-input", min);
    await utils.clearInputAndTypeInto("clip-max-input", max);
    await utils.clickOn("clip-commit");
  },

  async createCategory(categoryName) {
    await utils.clickOn("open-annotation-dialog");
    await utils.typeInto("new-category-name", categoryName);
    await utils.clickOn("submit-category");
  },

  async renameCategory(oldCatgoryName, newCategoryName) {
    await utils.clickOn(`${oldCatgoryName}:see-actions`);
    await utils.clickOn(`${oldCatgoryName}:edit-category-mode`);
    await utils.clearInputAndTypeInto(
      `${oldCatgoryName}:edit-category-name-text`,
      newCategoryName
    );
    await utils.clickOn(`${oldCatgoryName}:submit-category-edit`);
  },

  async deleteCategory(categoryName) {
    await utils.clickOn(`${categoryName}:see-actions`);
    await utils.clickOn(`${categoryName}:delete-category`);
  },

  async createLabel(categoryName, labelName) {
    await utils.clickOn(`${categoryName}:see-actions`);
    await utils.clickOn(`${categoryName}:add-new-label-to-category`);
    await utils.typeInto(`${categoryName}:new-label-name`, labelName);
    await utils.clickOn(`${categoryName}:submit-label`);
  },

  async deleteLabel(categoryName, labelName) {
    await this.expandCategory(categoryName);
    await utils.clickOn(`${categoryName}:${labelName}:see-actions`);
    await utils.clickOn(`${categoryName}:${labelName}:delete-label`);
  },

  async renameLabel(categoryName, oldLabelName, newLabelName) {
    await this.expandCategory(categoryName);
    await utils.clickOn(`${categoryName}:${oldLabelName}:see-actions`);
    await utils.clickOn(`${categoryName}:${oldLabelName}:edit-label`);
    await utils.clearInputAndTypeInto(
      `${categoryName}:${oldLabelName}:edit-label-name`,
      newLabelName
    );
    await utils.clickOn(`${categoryName}:${oldLabelName}:submit-label-edit`);
  },

  async addGeneToSearch(geneName) {
    await utils.typeInto("gene-search", geneName);
    await page.keyboard.press("Enter");
    await page.waitForSelector(`[data-testid='histogram-${geneName}']`);
  },

  async subset(coordinatesAsPercent) {
    // In order to deselect the selection after the subset, make sure we have some clear part
    // of the scatterplot we can click on
    assert(coordinatesAsPercent.x2 < 0.99 || coordinatesAsPercent.y2 < 0.99);
    const lassoSelection = await this.calcDragCoordinates(
      "layout-graph",
      coordinatesAsPercent
    );
    await this.drag(
      "layout-graph",
      lassoSelection.start,
      lassoSelection.end,
      true
    );
    await utils.clickOn("subset-button");
    const clearCoordinate = await this.calcCoordinate(
      "layout-graph",
      0.5,
      0.99
    );
    await this.clickOnCoordinate("layout-graph", clearCoordinate);
  },

  async setSellSet(cellSet, cellSetNum) {
    for (const selection of cellSet.filter(
      (sel) => sel.kind === "categorical"
    )) {
      await this.selectCategory(selection.metadata, selection.values, true);
    }
    await this.cellSet(cellSetNum);
  },

  async runDiffExp(cellSet1, cellSet2) {
    await this.setSellSet(cellSet1, 1);
    await this.setSellSet(cellSet2, 2);
    await utils.clickOn("diffexp-button");
  },

  async bulkAddGenes(geneNames) {
    await utils.clickOn("section-bulk-add");
    await utils.typeInto("input-bulk-add", geneNames.join(","));
    await page.keyboard.press("Enter");
  },
});
