import puppeteer from "puppeteer";

const appUrlBase = "http://localhost:3000";

let browser;
let page;
const browserViewport = { width: 1280, height: 960 };
jest.setTimeout(100000);

beforeAll(async () => {
  browser = await puppeteer.launch({
    headless: false,
    slowMo: 200,
    devtools: true
  });
  page = await browser.newPage();
  page.setViewport(browserViewport);
  page.on("console", msg => console.log("PAGE LOG:", msg.text()));
});

afterAll(() => {
  if (!process.env.DEBUG) {
    browser.close();
  }
});

const getOneElementInnerHTML = async function(selector) {
  let text = await page.$eval(selector, el => el.innerHTML);
  return text;
};

const drag = async function(el_box, start, end, lasso = false) {
  const x1 = el_box.content[0].x + start.x;
  const x2 = el_box.content[0].x + end.x;
  const y1 = el_box.content[0].y + start.y;
  const y2 = el_box.content[0].y + end.y;
  await page.mouse.move(x1, y1);
  await page.mouse.down();
  if (lasso) {
    await page.mouse.move(x2, y1, { steps: 50 });
    await page.mouse.move(x2, y2, { steps: 50 });
    await page.mouse.move(x1, y2, { steps: 50 });
    await page.mouse.move(x1, y1, { steps: 50 });
  } else {
    await page.mouse.move(x2, y2);
  }
  await page.mouse.up();
};

const lasso = async function(el_box, start, end) {
  const x1 = el_box.content[0].x + start.x;
  const x2 = el_box.content[0].x + end.x;
  const y1 = el_box.content[0].y + start.y;
  const y2 = el_box.content[0].y + end.y;

  const mouse = page.mouse;
  await mouse.click(x1, y1);
  await mouse.click(x2, y1);
  await mouse.click(x2, y2);
  await mouse.click(x1, y2);
  await mouse.click(x1, y1);
};

describe("did launch", () => {
  test("page launched", async () => {
    await page.goto(appUrlBase);
    let el = await getOneElementInnerHTML("[data-testid='header']");
    expect(el).toBe("cellxgene: pbmc3k");
  });
});

describe("search for genes", () => {
  test("search for known gene and add to metadata", async () => {
    await page.goto(appUrlBase);
    await page.waitForSelector("[ data-testid='gene-search']");
    // blueprint's  typeahead is treating typing weird, clicking first solves this
    await page.click("[data-testid='gene-search']");
    await page.type("[data-testid='gene-search']", "ACD");
    await page.focus("[data-testid='gene-search']");
    await page.keyboard.press("Enter");
    await page.waitForSelector("[data-testid='histogram-ACD']");
  });
});

describe("select cells and diffexp", () => {
  test.only("selects cells from layout and adds to cell set 1", async () => {
    await page.goto(appUrlBase);
    const layout = await page.waitForSelector("[data-testid='layout']");
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
    await lasso(size, cellset1.start, cellset1.end);
    await page.click("[data-testid='cellset_button_1");
    let button = await getOneElementInnerHTML("[data-testid='cellset_button_1");
    expect(button).toMatch(/26 cells/);
  });
  test("selects cells from layout and adds to cell set 2", async () => {
    await page.goto(appUrlBase);
    const layout = await page.waitForSelector("[data-testid='layout']");
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
    await lasso(size, cellset2.start, cellset2.end);
    await page.click("[data-testid='cellset_button_2");
    let button = await getOneElementInnerHTML("[data-testid='cellset_button_2");
    expect(button).toMatch(/49 cells/);
  });
  test("selects cells, saves them and performs diffexp", async () => {
    await page.goto(appUrlBase);
    const layout = await page.waitForSelector("[data-testid='layout']");
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
    await lasso(size, cellset1.start, cellset1.end);
    await page.click("[data-testid='cellset_button_1");
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
    await lasso(size, cellset2.start, cellset2.end);
    await page.click("[data-testid='cellset_button_2");
    await page.click("[data-testid='diffexp_button");
    await page.waitForSelector("[data-testclass='histogram-diffexp']");
    const diffexps = await page.$$eval(
      "[data-testclass='histogram-diffexp']",
      divs => {
        return divs.map(div =>
          div.id.substring("histogram_".length, div.id.length)
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
      "[data-testid='histogram_n_genes_svg_brush'] > .overlay"
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
    await drag(hist_size, draghist.start, draghist.end);
  });
});
