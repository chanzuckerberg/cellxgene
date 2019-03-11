import puppeteer from "puppeteer";
import { appUrlBase, DEBUG, DEV } from "./const";
import { puppeteerUtils } from "./utils";

let browser, page, utils;
const browserViewport = { width: 1280, height: 960 };

if (DEBUG) jest.setTimeout(100000);

beforeAll(async () => {
  const browser_params = DEV
    ? { headless: false, slowMo: 50 }
    : DEBUG
    ? { headless: false, slowMo: 200, devtools: true }
    : {};
  browser = await puppeteer.launch(browser_params);
  page = await browser.newPage();
  page.setViewport(browserViewport);
  if (DEV || DEBUG)
    page.on("console", msg => console.log("PAGE LOG:", msg.text()));
  utils = puppeteerUtils(page);
});

afterAll(() => {
  if (!DEBUG) {
    browser.close();
  }
});

describe("did launch", () => {
  test("page launched", async () => {
    await page.goto(appUrlBase);
    let el = await utils.getOneElementInnerHTML("[data-testid='header']");
    expect(el).toBe("cellxgene: pbmc3k");
  });
});

describe("search for genes", () => {
  test("search for known gene and add to metadata", async () => {
    await page.goto(appUrlBase);
    await page.waitForSelector("[ data-testid='gene-search']");
    // blueprint's  typeahead is treating typing weird, clicking & waiting first solves this
    await page.click("[data-testid='gene-search']");
    await page.waitFor(200);
    await page.type("[data-testid='gene-search']", "ACD");
    await page.keyboard.press("Enter");
    await page.waitForSelector("[data-testid='histogram-ACD']");
  });
});

describe("select cells and diffexp", () => {
  test("selects cells from layout and adds to cell set 1", async () => {
    await page.goto(appUrlBase);
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
    await utils.drag(size, cellset1.start, cellset1.end, true);
    await page.click("[data-testid='cellset-button-1");
    let button = await utils.getOneElementInnerHTML(
      "[data-testid='cellset-button-1"
    );
    expect(button).toMatch(/26 cells/);
  });

  test("selects cells from layout and adds to cell set 2", async () => {
    await page.goto(appUrlBase);
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
    await utils.drag(size, cellset2.start, cellset2.end, true);
    await page.click("[data-testid='cellset-button-2");
    let button = await utils.getOneElementInnerHTML(
      "[data-testid='cellset-button-2"
    );
    expect(button).toMatch(/49 cells/);
  });

  test("selects cells, saves them and performs diffexp", async () => {
    await page.goto(appUrlBase);
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
    await utils.drag(size, cellset1.start, cellset1.end, true);
    await page.click("[data-testid='cellset-button-1");
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
    await utils.drag(size, cellset2.start, cellset2.end, true);
    await page.click("[data-testid='cellset-button-2");
    await page.click("[data-testid='diffexp-button");
    await page.waitForSelector("[data-testclass='histogram-diffexp']");
    const diffexps = await page.$$eval(
      "[data-testclass='histogram-diffexp']",
      divs => {
        return divs.map(div =>
          div.id.substring("histogram-".length, div.id.length)
        );
      }
    );
    expect(diffexps).toMatchObject([
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
    await page.goto(appUrlBase);
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
    await utils.drag(hist_size, draghist.start, draghist.end);
  });
});
