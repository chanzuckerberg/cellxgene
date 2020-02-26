/*
Smoke test suite that will be run in Travis CI

Tests included in this file are expected to be relatively stable and test core features
 */
import { appUrlBase, DEBUG, DATASET } from "./config";
import { setupTestBrowser } from "./puppeteerUtils";
import { datasets } from "./data";

let browser, page, utils, cxgActions, spy;
const browserViewport = { width: 1280, height: 960 };
let data = datasets[DATASET];

beforeAll(async () => {
    [browser, page, utils, cxgActions] = await setupTestBrowser(browserViewport);
});

beforeEach(async () => {
  await page.goto(appUrlBase);
});

afterAll(() => {
  if (!DEBUG) {
    browser.close();
  }
});

describe("did launch", () => {
  test("page launched", async () => {
    let el = await utils.getOneElementInnerHTML("[data-testid='header']");
    expect(el).toBe(data.title);
  });
});

describe("metadata loads", () => {
  test("categories and values from dataset appear", async () => {
    for (const label in data.categorical) {
      await utils.waitByID(`category-${label}`);
      const categoryName = await utils.getOneElementInnerText(
        `[data-testid="category-${label}"]`
      );
      expect(categoryName).toMatch(label);
      await utils.clickOn(`${label}:category-expand`);
      const categories = await cxgActions.getAllCategoriesAndCounts(label);
      expect(Object.keys(categories)).toMatchObject(
        Object.keys(data.categorical[label])
      );
      expect(Object.values(categories)).toMatchObject(
        Object.values(data.categorical[label])
      );
    }
  });

  test("continuous data appears", async () => {
    for (const label in data.continuous) {
      await utils.waitByID(`histogram-${label}`);
    }
  });
});

describe("cell selection", () => {
  test("selects all cells cellset 1", async () => {
    const cellCount = await cxgActions.cellSet(1);
    expect(cellCount).toBe(data.dataframe.nObs);
  });

  test("selects all cells cellset 2", async () => {
    const cellCount = await cxgActions.cellSet(2);
    expect(cellCount).toBe(data.dataframe.nObs);
  });

  test("selects cells via lasso", async () => {
    for (const cellset of data.cellsets.lasso) {
      const cellset1 = await cxgActions.calcDragCoordinates(
        "layout-graph",
        cellset["coordinates-as-percent"]
      );
      await cxgActions.drag("layout-graph", cellset1.start, cellset1.end, true);
      const cellCount = await cxgActions.cellSet(1);
      expect(cellCount).toBe(cellset.count);
    }
  });

  test("selects cells via categorical", async () => {
    for (const cellset of data.cellsets.categorical) {
      await utils.clickOn(`${cellset.metadata}:category-expand`);
      await utils.clickOn(`${cellset.metadata}:category-select`);
      for (const val of cellset.values) {
        await utils.clickOn(
          `categorical-value-select-${cellset.metadata}-${val}`
        );
      }
      const cellCount = await cxgActions.cellSet(1);
      expect(cellCount).toBe(cellset.count);
    }
  });

  test("selects cells via continuous", async () => {
    for (const cellset of data.cellsets.continuous) {
      const histBrushableAreaId = `histogram-${cellset.metadata}-plot-brushable-area`;
      const coords = await cxgActions.calcDragCoordinates(
        histBrushableAreaId,
        cellset["coordinates-as-percent"]
      );
      await cxgActions.drag(histBrushableAreaId, coords.start, coords.end);
      const cellCount = await cxgActions.cellSet(1);
      expect(cellCount).toBe(cellset.count);
    }
  });
});

describe("gene entry", () => {
  test("search for single gene", async () => {
    // blueprint's  typeahead is treating typing weird, clicking & waiting first solves this
    await utils.typeInto("gene-search", data.genes.search);
    await page.keyboard.press("Enter");
    await page.waitForSelector(
      `[data-testid='histogram-${data.genes.search}']`
    );
  });

  test("bulk add genes", async () => {
    await cxgActions.reset();
    const testGenes = data.genes.bulkadd;
    await utils.clickOn("section-bulk-add");
    await utils.typeInto("input-bulk-add", testGenes.join(","));
    await page.keyboard.press("Enter");

    const allHistograms = await cxgActions.getAllHistograms(
      "histogram-user-gene",
      testGenes
    );
    expect(allHistograms).toEqual(expect.arrayContaining(testGenes));
    expect(allHistograms.length).toEqual(testGenes.length);
  });
});

describe("diffexp", () => {
  test("selects cells, saves them and performs diffexp", async () => {
    for (const select of data.diffexp.cellset1) {
      if (select.kind === "categorical") {
        await cxgActions.selectCategory(select.metadata, select.values, true);
      }
    }
    await cxgActions.cellSet(1);
    for (const select of data.diffexp.cellset2) {
      if (select.kind === "categorical") {
        await cxgActions.selectCategory(select.metadata, select.values, true);
      }
    }
    await cxgActions.cellSet(2);
    await utils.clickOn("diffexp-button");
    const allHistograms = await cxgActions.getAllHistograms(
      "histogram-diffexp",
      data.diffexp["gene-results"]
    );
    expect(allHistograms).toEqual(
      expect.arrayContaining(data.diffexp["gene-results"])
    );
    expect(allHistograms.length).toEqual(data.diffexp["gene-results"].length);
  });
});

describe("subset/reset", () => {
  test("subset - cell count matches", async () => {
    for (const select of data.subset.cellset1) {
      if (select.kind === "categorical") {
        await cxgActions.selectCategory(select.metadata, select.values, true);
      }
    }
    await utils.clickOn("subset-button");
    for (const label in data.subset.categorical) {
      const categories = await cxgActions.getAllCategoriesAndCounts(label);
      expect(Object.keys(categories)).toMatchObject(
        Object.keys(data.subset.categorical[label])
      );
      expect(Object.values(categories)).toMatchObject(
        Object.values(data.subset.categorical[label])
      );
    }
  });

  test("reset after subset", async () => {
    for (const select of data.subset.cellset1) {
      if (select.kind === "categorical") {
        await cxgActions.selectCategory(select.metadata, select.values, true);
      }
    }
    await utils.clickOn("subset-button");
    for (const label in data.subset.categorical) {
      const categories = await cxgActions.getAllCategoriesAndCounts(label);
      expect(Object.keys(categories)).toMatchObject(
        Object.keys(data.subset.categorical[label])
      );
      expect(Object.values(categories)).toMatchObject(
        Object.values(data.subset.categorical[label])
      );
    }
    await cxgActions.reset();
    for (const label in data.categorical) {
      await utils.waitByID(`category-${label}`);
      const categoryName = await utils.getOneElementInnerText(
        `[data-testid="category-${label}"]`
      );
      expect(categoryName).toMatch(label);
      const categories = await cxgActions.getAllCategoriesAndCounts(label);
      expect(Object.keys(categories)).toMatchObject(
        Object.keys(data.categorical[label])
      );
      expect(Object.values(categories)).toMatchObject(
        Object.values(data.categorical[label])
      );
    }
  });

  test("lasso after subset", async () => {
    for (const select of data.subset.cellset1) {
      if (select.kind === "categorical") {
        await cxgActions.selectCategory(select.metadata, select.values, true);
      }
    }
    await utils.clickOn("subset-button");
    const lassoSelection = await cxgActions.calcDragCoordinates(
      "layout-graph",
      data.subset.lasso["coordinates-as-percent"]
    );
    await cxgActions.drag(
      "layout-graph",
      lassoSelection.start,
      lassoSelection.end,
      true
    );
    const cellCount = await cxgActions.cellSet(1);
    expect(cellCount).toBe(data.subset.lasso.count);
  });
});

describe("scatter plot", () => {
  test("scatter plot appears", async () => {
    await cxgActions.reset();
    const testGenes = data.scatter.genes;
    await utils.clickOn("section-bulk-add");
    await utils.typeInto("input-bulk-add", Object.values(testGenes).join(","));
    await page.keyboard.press("Enter");
    await utils.clickOn(`plot-x-${data.scatter.genes.x}`);
    await utils.clickOn(`plot-y-${data.scatter.genes.y}`);
    await utils.waitByID("scatterplot");
  });
});

describe("clipping", () => {
  test("clip continuous", async () => {
    await cxgActions.clip(data.clip.min, data.clip.max);
    const histBrushableAreaId = `histogram-${data.clip.metadata}-plot-brushable-area`;
    const coords = await cxgActions.calcDragCoordinates(
      histBrushableAreaId,
      data.clip["coordinates-as-percent"]
    );
    await cxgActions.drag(histBrushableAreaId, coords.start, coords.end);
    const cellCount = await cxgActions.cellSet(1);
    expect(cellCount).toBe(data.clip.count);
  });

  test("clip gene", async () => {
    await utils.typeInto("gene-search", data.clip.gene);
    await page.keyboard.press("Enter");
    await page.waitForSelector(`[data-testid='histogram-${data.clip.gene}']`);
    await cxgActions.clip(data.clip.min, data.clip.max);
    const histBrushableAreaId = `histogram-${data.clip.gene}-plot-brushable-area`;
    const coords = await cxgActions.calcDragCoordinates(
      histBrushableAreaId,
      data.clip["coordinates-as-percent"]
    );
    await cxgActions.drag(histBrushableAreaId, coords.start, coords.end);
    const cellCount = await cxgActions.cellSet(1);
    expect(cellCount).toBe(data.clip["gene-cell-count"]);
  });
});

// interact with UI elements just that they do not break
describe("ui elements don't error", () => {
  test("color by", async () => {
    for (const label in data.categorical) {
      await utils.clickOn(`colorby-${label}`);
    }
    for (const label in data.continuous) {
      await utils.clickOn(`colorby-${label}`);
    }
  });

  test("color by for gene", async () => {
    await utils.typeInto("gene-search", data.genes.search);
    await page.keyboard.press("Enter");
    await page.waitForSelector(
      `[data-testid='histogram-${data.genes.search}']`
    );
    await utils.clickOn(`colorby-${data.genes.search}`);
  });

  test("pan and zoom", async () => {
    await utils.clickOn("mode-pan-zoom");
    const panCoords = await cxgActions.calcDragCoordinates(
      "layout-graph",
      data.pan["coordinates-as-percent"]
    );
    await cxgActions.drag(
      "layout-graph",
      panCoords.start,
      panCoords.end,
      false
    );
    await page.evaluate("window.scrollBy(0, 1000);");
  });
});
