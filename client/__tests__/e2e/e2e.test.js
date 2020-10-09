/**
 * Smoke test suite that will be run in Travis CI
 * Tests included in this file are expected to be relatively stable and test core features
 */

/* eslint-disable no-await-in-loop -- await in loop is needed to emulate sequential user actions  */
import { appUrlBase, DATASET } from "./config";

import { datasets } from "./data";

import {
  clickOn,
  getAllByClass,
  getElementCoordinates,
  getOneElementInnerHTML,
  getTestId,
  goToPage,
  typeInto,
  waitByID,
  clickOnUntil,
} from "./puppeteerUtils";

import {
  addGeneToSearch,
  bulkAddGenes,
  calcDragCoordinates,
  clip,
  drag,
  getAllCategoriesAndCounts,
  getAllHistograms,
  getCellSetCount,
  runDiffExp,
  selectCategory,
  subset,
  login,
  logout,
} from "./cellxgeneActions";

const data = datasets[DATASET];

describe("did launch", () => {
  test("page launched", async () => {
    await goToPage(appUrlBase);

    const element = await getOneElementInnerHTML(getTestId("header"));

    expect(element).toMatchSnapshot();
  });
});

describe("metadata loads", () => {
  test("categories and values from dataset appear", async () => {
    await goToPage(appUrlBase);

    for (const label of Object.keys(data.categorical)) {
      const element = await getOneElementInnerHTML(
        getTestId(`category-${label}`)
      );

      expect(element).toMatchSnapshot();

      await clickOn(`${label}:category-expand`);

      const categories = await getAllCategoriesAndCounts(label);

      expect(Object.keys(categories)).toMatchObject(
        Object.keys(data.categorical[label])
      );

      expect(Object.values(categories)).toMatchObject(
        Object.values(data.categorical[label])
      );
    }
  });

  test("continuous data appears", async () => {
    await goToPage(appUrlBase);

    for (const label of Object.keys(data.continuous)) {
      await waitByID(`histogram-${label}`);
    }
  });
});

describe("cell selection", () => {
  test("selects all cells cellset 1", async () => {
    await goToPage(appUrlBase);

    const cellCount = await getCellSetCount(1);
    expect(cellCount).toBe(data.dataframe.nObs);
  });

  test("selects all cells cellset 2", async () => {
    await goToPage(appUrlBase);

    const cellCount = await getCellSetCount(2);
    expect(cellCount).toBe(data.dataframe.nObs);
  });

  test("selects cells via lasso", async () => {
    await goToPage(appUrlBase);

    for (const cellset of data.cellsets.lasso) {
      const cellset1 = await calcDragCoordinates(
        "layout-graph",
        cellset["coordinates-as-percent"]
      );

      await drag("layout-graph", cellset1.start, cellset1.end, true);
      const cellCount = await getCellSetCount(1);
      expect(cellCount).toBe(cellset.count);
    }
  });

  test("selects cells via categorical", async () => {
    await goToPage(appUrlBase);

    for (const cellset of data.cellsets.categorical) {
      await clickOn(`${cellset.metadata}:category-expand`);
      await clickOn(`${cellset.metadata}:category-select`);

      for (const value of cellset.values) {
        await clickOn(`categorical-value-select-${cellset.metadata}-${value}`);
      }

      const cellCount = await getCellSetCount(1);

      expect(cellCount).toBe(cellset.count);
    }
  });

  test("selects cells via continuous", async () => {
    await goToPage(appUrlBase);

    for (const cellset of data.cellsets.continuous) {
      const histBrushableAreaId = `histogram-${cellset.metadata}-plot-brushable-area`;

      const coords = await calcDragCoordinates(
        histBrushableAreaId,
        cellset["coordinates-as-percent"]
      );

      await drag(histBrushableAreaId, coords.start, coords.end);

      const cellCount = await getCellSetCount(1);

      expect(cellCount).toBe(cellset.count);
    }
  });
});

describe("gene entry", () => {
  test("search for single gene", async () => {
    await goToPage(appUrlBase);
    await addGeneToSearch(data.genes.search);
  });

  test("bulk add genes", async () => {
    await goToPage(appUrlBase);

    const testGenes = data.genes.bulkadd;

    await bulkAddGenes(testGenes);
    const allHistograms = await getAllHistograms(
      "histogram-user-gene",
      testGenes
    );

    expect(allHistograms).toEqual(expect.arrayContaining(testGenes));
    expect(allHistograms).toHaveLength(testGenes.length);
  });
});

describe("differential expression", () => {
  test("selects cells, saves them and performs diffexp", async () => {
    await goToPage(appUrlBase);

    await runDiffExp(data.diffexp.cellset1, data.diffexp.cellset2);

    const allHistograms = await getAllHistograms(
      "histogram-diffexp",
      data.diffexp["gene-results"]
    );

    expect(allHistograms).toEqual(
      expect.arrayContaining(data.diffexp["gene-results"])
    );

    expect(allHistograms).toHaveLength(data.diffexp["gene-results"].length);
  });
});

describe("subset", () => {
  test("subset - cell count matches", async () => {
    await goToPage(appUrlBase);

    for (const select of data.subset.cellset1) {
      if (select.kind === "categorical") {
        await selectCategory(select.metadata, select.values, true);
      }
    }

    await clickOn("subset-button");

    for (const label of Object.keys(data.subset.categorical)) {
      const categories = await getAllCategoriesAndCounts(label);

      expect(Object.keys(categories)).toMatchObject(
        Object.keys(data.subset.categorical[label])
      );

      expect(Object.values(categories)).toMatchObject(
        Object.values(data.subset.categorical[label])
      );
    }
  });

  test("lasso after subset", async () => {
    await goToPage(appUrlBase);

    for (const select of data.subset.cellset1) {
      if (select.kind === "categorical") {
        await selectCategory(select.metadata, select.values, true);
      }
    }

    await clickOn("subset-button");

    const lassoSelection = await calcDragCoordinates(
      "layout-graph",
      data.subset.lasso["coordinates-as-percent"]
    );

    await drag("layout-graph", lassoSelection.start, lassoSelection.end, true);

    const cellCount = await getCellSetCount(1);
    expect(cellCount).toBe(data.subset.lasso.count);
  });

  test("undo selection appends the top diff exp genes to user defined genes", async () => {
    await goToPage(appUrlBase);

    const userDefinedGenes = data.genes.bulkadd;
    const diffExpGenes = data.diffexp["gene-results"];

    await bulkAddGenes(userDefinedGenes);
    const userDefinedHistograms = await getAllHistograms(
      "histogram-user-gene",
      userDefinedGenes
    );

    expect(userDefinedHistograms).toEqual(
      expect.arrayContaining(userDefinedGenes)
    );

    await subset({ x1: 0.15, y1: 0.1, x2: 0.98, y2: 0.98 });
    await runDiffExp(data.diffexp.cellset1, data.diffexp.cellset2);

    const diffExpHistograms = await getAllHistograms(
      "histogram-diffexp",
      diffExpGenes
    );

    expect(diffExpHistograms).toEqual(expect.arrayContaining(diffExpGenes));

    await clickOn("reset-subset-button");
    const expected = [].concat(userDefinedGenes, diffExpGenes);
    const userDefinedHistogramsAfterSubset = await getAllHistograms(
      "histogram-user-gene",
      expected
    );

    expect(userDefinedHistogramsAfterSubset).toEqual(
      expect.arrayContaining(expected)
    );
  });

  test("subset selection appends the top diff exp genes to user defined genes", async () => {
    await goToPage(appUrlBase);

    const userDefinedGenes = data.genes.bulkadd;
    const diffExpGenes = data.diffexp["gene-results"];

    await bulkAddGenes(userDefinedGenes);
    const userDefinedHistograms = await getAllHistograms(
      "histogram-user-gene",
      userDefinedGenes
    );

    expect(userDefinedHistograms).toEqual(
      expect.arrayContaining(userDefinedGenes)
    );

    await subset({ x1: 0.15, y1: 0.1, x2: 0.98, y2: 0.98 });
    await runDiffExp(data.diffexp.cellset1, data.diffexp.cellset2);

    const diffExpHistograms = await getAllHistograms(
      "histogram-diffexp",
      diffExpGenes
    );

    expect(diffExpHistograms).toEqual(expect.arrayContaining(diffExpGenes));

    await subset({ x1: 0.16, y1: 0.11, x2: 0.97, y2: 0.97 });

    const expected = [].concat(userDefinedGenes, diffExpGenes);
    const userDefinedHistogramsAfterSubset = await getAllHistograms(
      "histogram-user-gene",
      expected
    );

    expect(userDefinedHistogramsAfterSubset).toEqual(
      expect.arrayContaining(expected)
    );
  });
});

describe("scatter plot", () => {
  test("scatter plot appears", async () => {
    await goToPage(appUrlBase);

    await bulkAddGenes(Object.values(data.scatter.genes));
    await clickOn(`plot-x-${data.scatter.genes.x}`);
    await clickOn(`plot-y-${data.scatter.genes.y}`);
    await waitByID("scatterplot");
  });
});

describe("clipping", () => {
  test("clip continuous", async () => {
    await goToPage(appUrlBase);

    await clip(data.clip.min, data.clip.max);
    const histBrushableAreaId = `histogram-${data.clip.metadata}-plot-brushable-area`;
    const coords = await calcDragCoordinates(
      histBrushableAreaId,
      data.clip["coordinates-as-percent"]
    );
    await drag(histBrushableAreaId, coords.start, coords.end);
    const cellCount = await getCellSetCount(1);
    expect(cellCount).toBe(data.clip.count);
  });

  test("clip gene", async () => {
    await goToPage(appUrlBase);

    await typeInto("gene-search", data.clip.gene);
    await page.keyboard.press("Enter");

    await page.waitForSelector(`[data-testid='histogram-${data.clip.gene}']`);

    await clip(data.clip.min, data.clip.max);

    const histBrushableAreaId = `histogram-${data.clip.gene}-plot-brushable-area`;

    const coords = await calcDragCoordinates(
      histBrushableAreaId,
      data.clip["coordinates-as-percent"]
    );

    await drag(histBrushableAreaId, coords.start, coords.end);

    const cellCount = await getCellSetCount(1);

    expect(cellCount).toBe(data.clip["gene-cell-count"]);
  });
});

// interact with UI elements just that they do not break
describe("ui elements don't error", () => {
  test("color by", async () => {
    await goToPage(appUrlBase);

    const allLabels = [
      ...Object.keys(data.categorical),
      ...Object.keys(data.continuous),
    ];

    for (const label of allLabels) {
      await clickOn(`colorby-${label}`);
    }
  });

  test("color by for gene", async () => {
    await goToPage(appUrlBase);

    await typeInto("gene-search", data.genes.search);
    await page.keyboard.press("Enter");
    await page.waitForSelector(
      `[data-testid='histogram-${data.genes.search}']`
    );
    await clickOn(`colorby-${data.genes.search}`);
  });

  test("pan and zoom", async () => {
    await goToPage(appUrlBase);

    await clickOn("mode-pan-zoom");
    const panCoords = await calcDragCoordinates(
      "layout-graph",
      data.pan["coordinates-as-percent"]
    );

    await drag("layout-graph", panCoords.start, panCoords.end, false);

    await page.evaluate("window.scrollBy(0, 1000);");
  });
});

describe("centroid labels", () => {
  test("labels are created", async () => {
    await goToPage(appUrlBase);

    const labels = Object.keys(data.categorical);

    await clickOn(`colorby-${labels[0]}`);
    await clickOn("centroid-label-toggle");

    // Toggle colorby for each category and check to see if labels are generated
    for (let i = 0, { length } = labels; i < length; i += 1) {
      const label = labels[i];
      // first label is already enabled
      if (i !== 0) await clickOn(`colorby-${label}`);
      const generatedLabels = await getAllByClass("centroid-label");
      // Number of labels generated should be equal to size of the object
      expect(generatedLabels).toHaveLength(
        Object.keys(data.categorical[label]).length
      );
    }
  });
});

describe("graph overlay", () => {
  test("transform centroids correctly", async () => {
    await goToPage(appUrlBase);

    const category = Object.keys(data.categorical)[0];

    await clickOn(`colorby-${category}`);
    await clickOn("centroid-label-toggle");
    await clickOn("mode-pan-zoom");

    const panCoords = await calcDragCoordinates(
      "layout-graph",
      data.pan["coordinates-as-percent"]
    );

    const categoryValue = Object.keys(data.categorical[category])[0];
    const initialCoordinates = await getElementCoordinates(
      `${categoryValue}-centroid-label`
    );

    await drag("layout-graph", panCoords.start, panCoords.end, false);
    const terminalCoordinates = await getElementCoordinates(
      `${categoryValue}-centroid-label`
    );

    expect(terminalCoordinates[0] - initialCoordinates[0]).toBeCloseTo(
      panCoords.end.x - panCoords.start.x
    );
    expect(terminalCoordinates[1] - initialCoordinates[1]).toBeCloseTo(
      panCoords.end.y - panCoords.start.y
    );
  });
});

test("pan zoom mode resets lasso selection", async () => {
  await goToPage(appUrlBase);

  const panzoomLasso = data.features.panzoom.lasso;

  const lassoSelection = await calcDragCoordinates(
    "layout-graph",
    panzoomLasso["coordinates-as-percent"]
  );

  await drag("layout-graph", lassoSelection.start, lassoSelection.end, true);
  await waitByID("lasso-element", { visible: true });

  const initialCount = await getCellSetCount(1);

  expect(initialCount).toBe(panzoomLasso.count);

  await clickOn("mode-pan-zoom");
  await clickOn("mode-lasso");

  const modeSwitchCount = await getCellSetCount(1);

  expect(modeSwitchCount).toBe(initialCount);
});

test("lasso moves after pan", async () => {
  await goToPage(appUrlBase);

  const panzoomLasso = data.features.panzoom.lasso;
  const coordinatesAsPercent = panzoomLasso["coordinates-as-percent"];

  const lassoSelection = await calcDragCoordinates(
    "layout-graph",
    coordinatesAsPercent
  );

  await drag("layout-graph", lassoSelection.start, lassoSelection.end, true);
  await waitByID("lasso-element", { visible: true });

  const initialCount = await getCellSetCount(1);

  expect(initialCount).toBe(panzoomLasso.count);

  await clickOn("mode-pan-zoom");

  const panCoords = await calcDragCoordinates(
    "layout-graph",
    coordinatesAsPercent
  );

  await drag("layout-graph", panCoords.start, panCoords.end, false);
  await clickOn("mode-lasso");

  const panCount = await getCellSetCount(2);

  expect(panCount).toBe(initialCount);
});

const describeIfCalledByMakeFileTarget =
  process.env.CXG_AUTH_TYPE === "test" ? describe : describe.skip;

describeIfCalledByMakeFileTarget("auth buttons", () => {
  test("login then logout", async () => {
    await goToPage(appUrlBase);
    await clickOnUntil("log-in", async () => {
      await page.waitForNavigation({ waitUntil: "networkidle0" });
      await waitByID("user-info");
    });
    await logout();
  });
});

const conditionalDescribe =
  process.env.TEST_AUTH_INTEGRATION === "true" ? describe : describe.skip;

conditionalDescribe("AuthN Integration", () => {
  it("logs in", async () => {
    await login();
  });

  it("logs out", async () => {
    await login();
    await logout();
  });
});
/* eslint-enable no-await-in-loop -- await in loop is needed to emulate sequential user actions */
