import puppeteer from "puppeteer";
import { appUrlBase, DEBUG, DEV } from "./const";
import { puppeteerUtils, cellxgeneActions } from "./utils";

let browser, page, utils, cxgActions;
const browserViewport = { width: 1280, height: 960 };

if (DEBUG) jest.setTimeout(100000);
if (DEV) jest.setTimeout(10000);

// TODO are tests robust?

beforeAll(async () => {
  const browser_params = DEV
    ? { headless: false, slowMo: 10 }
    : DEBUG
    ? { headless: false, slowMo: 200, devtools: true }
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

describe("did launch", () => {
  test("page launched", async () => {
    let el = await utils.getOneElementInnerHTML("[data-testid='header']");
    expect(el).toBe("cellxgene: pbmc3k");
  });
});

describe("search for genes", () => {
  test("search for known gene and add to metadata", async () => {
    // blueprint's  typeahead is treating typing weird, clicking & waiting first solves this
    await utils.typeInto("gene-search", "ACD");
    await page.keyboard.press("Enter");
    await page.waitForSelector("[data-testid='histogram-ACD']");
  });
});

describe("select cells and diffexp", () => {
  test("selects cells from layout and adds to cell set 1", async () => {
    const layout = await page.waitForSelector("[data-testid='layout-graph']");
    const size = await layout.boxModel();
    const cellset1 = {
      start: {
        x: Math.floor(size.width * 0.25),
        y: Math.floor(size.height * 0.25)
      },
      end: {
        x: Math.floor(size.width * 0.35),
        y: Math.floor(size.height * 0.35)
      }
    };
    await cxgActions.drag(size, cellset1.start, cellset1.end, true);
    const cell_count = await cxgActions.cellSet(1)
    expect(cell_count).toBe("26")
  });

  test("selects cells from layout and adds to cell set 2", async () => {
    const layout = await page.waitForSelector("[data-testid='layout-graph']");
    const size = await layout.boxModel();
    const cellset2 = {
      start: {
        x: Math.floor(size.width * 0.45),
        y: Math.floor(size.height * 0.45)
      },
      end: {
        x: Math.floor(size.width * 0.55),
        y: Math.floor(size.height * 0.55)
      }
    };
    await cxgActions.drag(size, cellset2.start, cellset2.end, true);
    const cell_count = await cxgActions.cellSet(2)
    expect(cell_count).toBe("49")
  });

  test("selects cells, saves them and performs diffexp", async () => {
    const layout = await page.waitForSelector("[data-testid='layout-graph']");
    const size = await layout.boxModel();
    const cellset1 = {
      start: {
        x: Math.floor(size.width * 0.25),
        y: Math.floor(size.height * 0.25)
      },
      end: {
        x: Math.floor(size.width * 0.35),
        y: Math.floor(size.height * 0.35)
      }
    };
    await cxgActions.drag(size, cellset1.start, cellset1.end, true);
    await utils.clickOn("cellset-button-1");
    const cellset2 = {
      start: {
        x: Math.floor(size.width * 0.45),
        y: Math.floor(size.height * 0.45)
      },
      end: {
        x: Math.floor(size.width * 0.55),
        y: Math.floor(size.height * 0.55)
      }
    };
    await cxgActions.drag(size, cellset2.start, cellset2.end, true);
    await utils.clickOn("cellset-button-2");
     await utils.clickOn("diffexp-button");
    const diffExpHists = await cxgActions.getAllHistograms("histogram-diffexp");
    expect(diffExpHists).toMatchObject([
      "HLA-DPA1",
      "HLA-DQA1",
      "HLA-DRB1",
      "HLA-DMA",
      "CST3",
      "HLA-DPB1",
      "HLA-DQB1",
      "LGALS2",
      "FCER1A",
      "LTB"
    ]);
  });
});

describe("brushable histogram", () => {
  test("can brush historgram", async () => {
    const hist = await page.waitForSelector(
      "[data-testid='histogram_n_genes_svg-brush'] > .overlay"
    );
    const hist_size = await hist.boxModel();
    const draghist = {
      start: {
        x: Math.floor(hist_size.width * 0.25),
        y: Math.floor(hist_size.height * 0.5)
      },
      end: {
        x: Math.floor(hist_size.width * 0.55),
        y: Math.floor(hist_size.height * 0.5)
      }
    };
    await cxgActions.drag(hist_size, draghist.start, draghist.end);
    const cell_count = await cxgActions.cellSet(1)
    expect(cell_count).toBe("1537")
  });
});

describe("bulk add genes", () => {
  test("add several genes in dataset and display them", async () => {
    await cxgActions.reset();
    const testGenes = ["S100A8", "FCGR3A", "LGALS2", "GSTP1"];
    await page.click("[data-testid='section-bulk-add']");
    await utils.typeInto("input-bulk-add", testGenes.join(","));
    await page.keyboard.press("Enter");
    const userGeneHist = await cxgActions.getAllHistograms(
      "histogram-user-gene"
    );
    expect(userGeneHist).toMatchObject(testGenes);
  });
});

describe("categorical data", () => {
  test("categories and values from dataset appear", async () => {
    await utils.waitByID("category-louvain");
    const louvain = await utils.getOneElementInnerText(
      '[data-testid="category-louvain"]'
    );
    expect(louvain).toMatch("louvain");
    await utils.clickOn("category-expand-louvain");
    const categories = await cxgActions.getAllCategoriesAndCounts("louvain");
    expect(Object.keys(categories)).toMatchObject([
      "B cells",
      "CD14+ Monocytes",
      "CD4 T cells",
      "CD8 T cells",
      "Dendritic cells",
      "FCGR3A+ Monocytes",
      "Megakaryocytes",
      "NK cells"
    ]);
    expect(Object.values(categories)).toMatchObject([
      "342",
      "480",
      "1144",
      "316",
      "37",
      "150",
      "15",
      "154"
    ]);
  })

  test("cell selection by categorical metadata", async () => {
    await utils.clickOn("category-expand-louvain");
    await utils.clickOn("category-select-louvain");
    await utils.clickOn("categorical-value-select-B cells");
    await utils.clickOn("categorical-value-select-Megakaryocytes");
    const cell_count = await cxgActions.cellSet(1)
    expect(cell_count).toBe("357")
  });
});
