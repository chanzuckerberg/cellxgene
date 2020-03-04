import {DEBUG, DEV} from "./config";
import puppeteer from "puppeteer";
import { strict as assert } from "assert";

export const puppeteerUtils = puppeteerPage => ({
  async waitByID(testId, props = {}) {
    return await puppeteerPage.waitForSelector(
      `[data-testid='${testId}']`,
      props
    );
  },

  async waitByClass(testClass, props = {}) {
    return await puppeteerPage.waitForSelector(
      `[data-testclass='${testClass}']`,
      props
    );
  },

  async waitForAllByIds(testIds) {
    await Promise.all(
      testIds.map(testId => puppeteerPage.waitForSelector(`[data-testid='${testId}']`))
    );
  },

  async getAllByClass(testClass) {
    return puppeteerPage.$$eval(
      `[data-testclass=${testClass}]`,
      eles => eles.map(ele => ele.dataset.testid)
    );
  },

  async typeInto(testid, text) {
    // blueprint's  typeahead is treating typing weird, clicking & waiting first solves this
    // only works for text without special characters
    await this.waitByID(testid);
    const selector = `[data-testid='${testid}']`;
    // type ahead can be annoying if you don't pause before you type
    await puppeteerPage.click(selector);
    await puppeteerPage.waitFor(200);
    await puppeteerPage.type(selector, text);
  },

  async clearInputAndTypeInto(testid, text) {
    await this.waitByID(testid);
    const selector = `[data-testid='${testid}']`;
    // only works for text without special characters
    // type ahead can be annoying if you don't pause before you type
    await puppeteerPage.click(selector);
    await puppeteerPage.waitFor(200);
    // select all

    await puppeteerPage.click(selector, { clickCount: 3 });
    await puppeteerPage.keyboard.press("Backspace");
    await puppeteerPage.type(selector, text);
  },

  async clickOn(testid, options={}) {
    await this.waitByID(testid);
    await puppeteerPage.click(`[data-testid='${testid}']`, options);
    await puppeteerPage.waitFor(50);
  },

  async getOneElementInnerHTML(selector) {
    await puppeteerPage.waitForSelector(selector);
    let text = await puppeteerPage.$eval(selector, el => el.innerHTML);
    return text;
  },

  async getOneElementInnerText(selector) {
    await puppeteerPage.waitForSelector(selector);
    let text = await puppeteerPage.$eval(selector, el => el.innerText);
    return text;
  }
});

export const cellxgeneActions = puppeteerPage => ({

  async drag(testid, start, end, lasso = false) {
    const layout = await puppeteerUtils(puppeteerPage).waitByID(testid);
    const elBox = await layout.boxModel();
    const x1 = elBox.content[0].x + start.x;
    const x2 = elBox.content[0].x + end.x;
    const y1 = elBox.content[0].y + start.y;
    const y2 = elBox.content[0].y + end.y;
    await puppeteerPage.mouse.move(x1, y1);
    await puppeteerPage.mouse.down();
    if (lasso) {
      await puppeteerPage.mouse.move(x2, y1);
      await puppeteerPage.mouse.move(x2, y2);
      await puppeteerPage.mouse.move(x1, y2);
      await puppeteerPage.mouse.move(x1, y1);
    } else {
      await puppeteerPage.mouse.move(x2, y2);
    }
    await puppeteerPage.mouse.up();
  },

  async clickOnCoordinate(testid, coord) {
    const layout = await puppeteerUtils(puppeteerPage).waitByID(testid);
    const elBox = await layout.boxModel();
    const x = elBox.content[0].x + coord.x;
    const y = elBox.content[0].y + coord.y;
    await puppeteerPage.mouse.click(x, y);
  },

  async getAllHistograms(testclass, testids) {
    const histTestIds = testids.map(tid => `histogram-${tid}`);
    // these load asynchronously, so we need to wait for each histogram individually
    await puppeteerUtils(puppeteerPage).waitForAllByIds(histTestIds);
    const allHistograms = await puppeteerUtils(puppeteerPage).getAllByClass(
      testclass
    );
    return allHistograms.map(hist => hist.replace(/^histogram-/, ""));
  },

  async getAllCategoriesAndCounts(category) {
    await puppeteerUtils(puppeteerPage).waitByClass("categorical-row");
    const categories = await puppeteerPage.$$eval(
      `[data-testid="category-${category}"] [data-testclass='categorical-row']`,
      els => {
        let result = {};
        els.forEach(el => {
          const cat = el.querySelector("[data-testclass='categorical-value']")
            .innerText;
          const count = el.querySelector(
            "[data-testclass='categorical-value-count']"
          ).innerText;
          result[cat] = count;
        });
        return result;
      }
    );
    return categories;
  },

  async cellSet(num) {
    await puppeteerUtils(puppeteerPage).clickOn(`cellset-button-${num}`);
    return await puppeteerUtils(puppeteerPage).getOneElementInnerText(
      `[data-testid='cellset-count-${num}']`
    );
  },

  async resetCategory(category) {
    const checkboxId = `${category}:category-select`;
    const utils = puppeteerUtils(puppeteerPage);
    await utils.waitByID(checkboxId);
    const checkedPseudoclass = await puppeteerPage.$eval(
      `[data-testid='${checkboxId}']`,
      el => {
        return el.matches(":checked");
      }
    );
    if (!checkedPseudoclass) {
      await utils.clickOn(checkboxId);
    }
    try {
      const categoryRow = await utils.waitByID(
        `${category}:category-expand`
      );
      const isExpanded = await categoryRow.$(
        "[data-testclass='category-expand-is-expanded']"
      );
      if (isExpanded) {
        await utils.clickOn(
          `${category}:category-expand`
        );
      }
    } catch {}
  },

  async calcCoordinate(testid, xAsPercent, yAsPercent) {
    const el = await puppeteerUtils(puppeteerPage).waitByID(testid);
    const size = await el.boxModel();
    return {
      x: Math.floor(size.width * xAsPercent),
      y: Math.floor(size.height * yAsPercent)
    }
  },

  async calcDragCoordinates(testid, coordinateAsPercent) {
    const coords = {
      start: await this.calcCoordinate(testid, coordinateAsPercent.x1, coordinateAsPercent.y1),
      end: await this.calcCoordinate(testid, coordinateAsPercent.x2, coordinateAsPercent.y2)
    };
    return coords;
  },

  async selectCategory(category, values, reset = true) {
    if (reset) await this.resetCategory(category);
    const utils = puppeteerUtils(puppeteerPage);
    await utils.clickOn(`${category}:category-expand`);
    await utils.clickOn(`${category}:category-select`);
    for (const val of values) {
      await utils.clickOn(`categorical-value-select-${category}-${val}`);
    }
  },

  async expandCategory(category) {
    const utils = puppeteerUtils(puppeteerPage);
    const expand = await utils.waitByID(`${category}:category-expand`);
    const notExpanded = await expand.$("[data-testclass='category-expand-is-not-expanded']");
    if (notExpanded) {
      await utils.clickOn(`${category}:category-expand`);
    }
  },

  async clip(min = 0, max = 100) {
    const utils = puppeteerUtils(puppeteerPage);
    await utils.clickOn("visualization-settings");
    await utils.clearInputAndTypeInto(
      "clip-min-input",
      min
    );
    await utils.clearInputAndTypeInto(
      "clip-max-input",
      max
    );
    await utils.clickOn("clip-commit");
  },

  async createCategory(categoryName) {
    const utils = puppeteerUtils(puppeteerPage);
    await utils.clickOn("open-annotation-dialog");
    await utils.typeInto("new-category-name", categoryName);
    await utils.clickOn("submit-category");
  },

  async renameCategory(oldCatgoryName, newCategoryName) {
    const utils = puppeteerUtils(puppeteerPage);
    await utils.clickOn(`${oldCatgoryName}:see-actions`);
    await utils.clickOn(`${oldCatgoryName}:edit-category-mode`);
    await utils.clearInputAndTypeInto(`${oldCatgoryName}:edit-category-name-text`, newCategoryName);
    await utils.clickOn(`${oldCatgoryName}:submit-category-edit`);
  },

  async deleteCategory(categoryName) {
    const utils = puppeteerUtils(puppeteerPage);
    await utils.clickOn(`${categoryName}:see-actions`);
    await utils.clickOn(`${categoryName}:delete-category`);
  },

  async createLabel(categoryName, labelName) {
    const utils = puppeteerUtils(puppeteerPage);
    await utils.clickOn(`${categoryName}:see-actions`);
    await utils.clickOn(`${categoryName}:add-new-label-to-category`);
    await utils.typeInto(`${categoryName}:new-label-name`, labelName);
    await utils.clickOn(`${categoryName}:submit-label`);
  },

  async deleteLabel(categoryName, labelName) {
    await this.expandCategory(categoryName);
    const utils = puppeteerUtils(puppeteerPage);
    await utils.clickOn(`${categoryName}:${labelName}:see-actions`);
    await utils.clickOn( `${categoryName}:${labelName}:delete-label`);
  },

  async renameLabel(categoryName, oldLabelName, newLabelName) {
    await this.expandCategory(categoryName);
    const utils = puppeteerUtils(puppeteerPage);
    await utils.clickOn(`${categoryName}:${oldLabelName}:see-actions`);
    await utils.clickOn(`${categoryName}:${oldLabelName}:edit-label`);
    await utils.clearInputAndTypeInto(
      `${categoryName}:${oldLabelName}:edit-label-name`,
      newLabelName
    );
    await puppeteerUtils(puppeteerPage).clickOn(`${categoryName}:${oldLabelName}:submit-label-edit`);
  },

  async addGeneToSearch(geneName) {
    await puppeteerUtils(puppeteerPage).typeInto("gene-search", geneName);
    await puppeteerPage.keyboard.press("Enter");
    await puppeteerPage.waitForSelector(
      `[data-testid='histogram-${geneName}']`
    );
  },

  async subset(coordinatesAsPercent) {
    // In order to deselect the selection after the subset, make sure we have some clear part
    // of the scatterplot we can click on
    assert(coordinatesAsPercent.x2 < 0.99 || coordinatesAsPercent.y2 < 0.99);
    const lassoSelection = await this.calcDragCoordinates( "layout-graph", coordinatesAsPercent);
    await this.drag("layout-graph", lassoSelection.start, lassoSelection.end, true );
    await puppeteerUtils(puppeteerPage).clickOn("subset-button");
    const clearCoordinate = await this.calcCoordinate(
      "layout-graph",
      0.5,
      0.99
    );
    await this.clickOnCoordinate("layout-graph", clearCoordinate);
  },

  async setSellSet(cellSet, cellSetNum) {
    for (const selection of cellSet) {
      if (selection.kind === "categorical") {
        await this.selectCategory(selection.metadata, selection.values, true);
      }
    }
    await this.cellSet(cellSetNum);
  },

  async runDiffExp(cellSet1, cellSet2) {
    await this.setSellSet(cellSet1, 1);
    await this.setSellSet(cellSet2, 2);
    await puppeteerUtils(puppeteerPage).clickOn("diffexp-button");
  },

  async bulkAddGenes(geneNames) {
    const utils = puppeteerUtils(puppeteerPage);
    await utils.clickOn("section-bulk-add");
    await utils.typeInto("input-bulk-add", geneNames.join(","));
    await puppeteerPage.keyboard.press("Enter");
  }
});

export async function setupTestBrowser(browserViewport) {
  const browserParams = DEV
      ? { headless: false, slowMo: 5 }
      : DEBUG
          ? { headless: false, slowMo: 100, devtools: true }
          : {};
  const browser = await puppeteer.launch(browserParams);
  const page = await browser.newPage();
  await page.setViewport(browserViewport);
  if (DEV || DEBUG) {
    page.on("console", async msg => {
      // If there is a console.error but an error is not thrown, this will ensure the test fails
      if (msg.type() === "error") {
        const errorMsgText = await Promise.all(
            // TODO can we do this without internal properties?
            msg.args().map(arg => arg._remoteObject.description)
        );
        throw new Error(`Console error: ${errorMsgText}`);
      }
      console.log(`PAGE LOG: ${msg.text()}`);
    });
  }
  page.on("pageerror", err => {
    throw new Error(`Console error: ${err}`);
  });
  const utils = puppeteerUtils(page);
  const cxgActions = cellxgeneActions(page);
  return [browser, page, utils, cxgActions];
}
