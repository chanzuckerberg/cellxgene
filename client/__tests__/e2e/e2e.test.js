import puppeteer from "puppeteer";
import { appUrlBase, DEBUG, DEV, DATASET } from "./config";
import { puppeteerUtils, cellxgeneActions } from "./puppeteer-utils";
import { datasets } from "./data";

let browser, page, utils, cxgActions, spy;
const browserViewport = { width: 1280, height: 960 };
let data = datasets[DATASET];

if (DEBUG) jest.setTimeout(100000);
if (DEV) jest.setTimeout(10000);

beforeAll(async () => {
  const browser_params = DEV
    ? { headless: false, slowMo: 5 }
    : DEBUG
    ? { headless: false, slowMo: 100, devtools: true }
    : {};
  browser = await puppeteer.launch(browser_params);
  page = await browser.newPage();
  await page.setViewport(browserViewport);
  if (DEV || DEBUG) {
    page.on("console", msg => console.log(`PAGE LOG: ${msg.text()}`));
  }
  page.on("pageerror", err => {
    throw new Error(`Console error: ${err}`);
  });
  utils = puppeteerUtils(page);
  cxgActions = cellxgeneActions(page);
});

beforeEach(async () => {
  await page.goto(appUrlBase);
});

afterAll(() => {
  if (!DEBUG) {
    browser.close();
  }
});

describe("did launch", async () => {
  test("page launched", async () => {
    let el = await utils.getOneElementInnerHTML("[data-testid='header']");
    expect(el).toBe(data.title);
  });
});

describe("metadata loads", async () => {
  test("categories and values from dataset appear", async () => {
    for (const label in data.categorical) {
      await utils.waitByID(`category-${label}`);
      const category_name = await utils.getOneElementInnerText(
        `[data-testid="category-${label}"]`
      );
      expect(category_name).toMatch(label);
      await utils.clickOn(`category-expand-${label}`);
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

describe("cell selection", async () => {
  test("selects all cells cellset 1", async () => {
    const cell_count = await cxgActions.cellSet(1);
    expect(cell_count).toBe(data.dataframe.nObs);
  });

  test("selects all cells cellset 2", async () => {
    const cell_count = await cxgActions.cellSet(2);
    expect(cell_count).toBe(data.dataframe.nObs);
  });

  test("selects cells via lasso", async () => {
    for (const cellset of data.cellsets.lasso) {
      const cellset1 = await cxgActions.calcDragCoordinates(
        "layout-graph",
        cellset["coordinates-as-percent"]
      );
      await cxgActions.drag("layout-graph", cellset1.start, cellset1.end, true);
      const cell_count = await cxgActions.cellSet(1);
      expect(cell_count).toBe(cellset.count);
    }
  });

  test("selects cells via categorical", async () => {
    for (const cellset of data.cellsets.categorical) {
      await utils.clickOn(`category-expand-${cellset.metadata}`);
      await utils.clickOn(`category-select-${cellset.metadata}`);
      for (const val of cellset.values) {
        await utils.clickOn(
          `categorical-value-select-${cellset.metadata}-${val}`
        );
      }
      const cell_count = await cxgActions.cellSet(1);
      expect(cell_count).toBe(cellset.count);
    }
  });

  test("selects cells via continuous", async () => {
    for (const cellset of data.cellsets.continuous) {
      const hist_id = `histogram-${cellset.metadata}-plot-brush`;
      const coords = await cxgActions.calcDragCoordinates(
        hist_id,
        cellset["coordinates-as-percent"]
      );
      await cxgActions.drag(hist_id, coords.start, coords.end);
      const cell_count = await cxgActions.cellSet(1);
      expect(cell_count).toBe(cellset.count);
    }
  });
});

describe("gene entry", async () => {
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
    const userGeneHist = await cxgActions.getAllHistograms(
      "histogram-user-gene"
    );
    expect(userGeneHist).toEqual(expect.arrayContaining(testGenes));
  });
});

describe("diffexp", async () => {
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
    const diffExpHists = await cxgActions.getAllHistograms("histogram-diffexp");
    expect(diffExpHists).toEqual(
      expect.arrayContaining(data.diffexp["gene-results"])
    );
  });
});
//

describe("subset/reset", async () => {
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
      const category_name = await utils.getOneElementInnerText(
        `[data-testid="category-${label}"]`
      );
      expect(category_name).toMatch(label);
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
    const lasso_selection = await cxgActions.calcDragCoordinates(
      "layout-graph",
      data.subset.lasso["coordinates-as-percent"]
    );
    await cxgActions.drag(
      "layout-graph",
      lasso_selection.start,
      lasso_selection.end,
      true
    );
    const cell_count = await cxgActions.cellSet(1);
    expect(cell_count).toBe(data.subset.lasso.count);
  });
});

describe("scatter plot", async () => {
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

// interact with UI elements just that they do not break
describe("ui elements don't error", async () => {
  test("color by", async () => {
    for (const label in data.categorical) {
      await utils.clickOn(`colorby-${label}`);
    }
    for (const label in data.continuous) {
      await utils.clickOn(`colorby-${label}`);
    }
  });

  test("pan and zoom", async () => {
    await utils.clickOn("mode-pan-zoom");
    const pan_coords = await cxgActions.calcDragCoordinates(
      "layout-graph",
      data.pan["coordinates-as-percent"]
    );
    await cxgActions.drag(
      "layout-graph",
      pan_coords.start,
      pan_coords.end,
      false
    );
    await page.evaluate(`window.scrollBy(0, 1000);`);
  });
});
