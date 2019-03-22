import puppeteer from "puppeteer";
import { appUrlBase, DEBUG, DEV, DATASET } from "./config";
import { puppeteerUtils, cellxgeneActions } from "./puppeteer-utils";
import { datasets } from "./data";

let browser, page, utils, cxgActions, store;
const browserViewport = { width: 1280, height: 960 };
let data = datasets[DATASET];

if (DEBUG) jest.setTimeout(100000);
if (DEV) jest.setTimeout(10000);

// TODO are tests robust?

beforeAll(async () => {
  const browser_params = DEV
    ? { headless: false, slowMo: 5 }
    : DEBUG
    ? { headless: false, slowMo: 100, devtools: true }
    : {};
  browser = await puppeteer.launch(browser_params);
  page = await browser.newPage();
  page.setViewport(browserViewport);
  if (DEV || DEBUG)
    page.on("console", msg => console.log("PAGE LOG:", msg.text()));
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
    expect(cell_count).toBe("2638");
  });

  test("selects all cells cellset 2", async () => {
    const cell_count = await cxgActions.cellSet(2);
    expect(cell_count).toBe(data.dataframe.nObs);
  });

  test("selects cells via lasso", async () => {
    for (let i = 0; i < data.cellsets.lasso.length; i++) {
      const cellset = data.cellsets.lasso[i];
      const cellset1 = await cxgActions.calcDragCoordinates(
        "layout-graph",
        cellset["coordinates-as-percent"]
      );
      await cxgActions.drag("layout-graph", cellset1.start, cellset1.end, true);
      const cell_count = await cxgActions.cellSet(1);
      expect(cell_count).toBe(cellset["count"]);
    }
  });

  test("selects cells via categorical", async () => {
    for (let i = 0; i < data.cellsets.categorical.length; i++) {
      const cellset = data.cellsets.categorical[i];
      await utils.clickOn(`category-expand-${cellset.metadata}`);
      await utils.clickOn(`category-select-${cellset.metadata}`);
      for (let i = 0; i < cellset.values.length; i++) {
        const val = cellset.values[i];
        await utils.clickOn(
          `categorical-value-select-${cellset.metadata}-${val}`
        );
      }
      const cell_count = await cxgActions.cellSet(1);
      expect(cell_count).toBe(cellset.count);
    }
  });

  test("selects cells via continuous", async () => {
    for (let i = 0; i < data.cellsets.continuous.length; i++) {
      const cellset = data.cellsets.continuous[i];
      const hist_id = `histogram_${cellset.metadata}_svg-brush`;
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
    const testGenes = data.genes["bulk add"];
    await page.click("[data-testid='section-bulk-add']");
    await utils.typeInto("input-bulk-add", testGenes.join(","));
    await page.keyboard.press("Enter");
    const userGeneHist = await cxgActions.getAllHistograms(
      "histogram-user-gene"
    );
    expect(userGeneHist).toMatchObject(testGenes);
  });
});

describe("diffexp", async () => {
  test("selects cells, saves them and performs diffexp", async () => {
    for (let i = 0; i < data.diffexp.cellset1.length; i++) {
      const select = data.diffexp.cellset1[i];
      if (select.kind === "categorical") {
        await cxgActions.selectCategory(select.metadata, select.values, false);
      }
    }
    const cellSet1 = await cxgActions.cellSet(1);
    expect(cellSet1).toBe("342");
    for (let i = 0; i < data.diffexp.cellset2.length; i++) {
      const select = data.diffexp.cellset2[i];
      if (select.kind === "categorical") {
        await cxgActions.selectCategory(select.metadata, select.values, true);
      }
    }

    const cellSet2 = await cxgActions.cellSet(2);
    expect(cellSet2).toBe("1298");

    await utils.clickOn("diffexp-button");
    const diffExpHists = await cxgActions.getAllHistograms("histogram-diffexp");
    expect(diffExpHists).toMatchObject(data.diffexp["gene-results"]);
  });
});
//

// describe("subset/reset", () => {});
//
// describe("scatter plot", () => {});
//
// describe("ui elements screen", () => {});
//
// describe("store test", () => {
//   test("this is a test", async () => {
//     store = await page.evaluate(() => {
//       window.cxg_store.getState();
//     });
//     console.log("Store2:", store);
//   });
// });
//
