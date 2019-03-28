export const puppeteerUtils = puppeteerPage => ({
  async waitByID(testid, props = {}) {
    return await puppeteerPage.waitForSelector(
      `[data-testid='${testid}']`,
      props
    );
  },

  async waitByClass(testclass, props = {}) {
    return await puppeteerPage.waitForSelector(
      `[data-testclass='${testclass}']`,
      props
    );
  },

  async typeInto(testid, text) {
    // only works for text without special characters
    await this.waitByID(testid);
    // type ahead can be annoying if you don't pause before you type
    await puppeteerPage.click(`[data-testid='${testid}']`);
    await puppeteerPage.waitFor(200);
    await puppeteerPage.type(`[data-testid='${testid}']`, text);
  },

  async clickOn(testid) {
    await this.waitByID(testid);
    await puppeteerPage.click(`[data-testid='${testid}']`);
    await puppeteerPage.waitFor(50);
  },

  async getOneElementInnerHTML(selector) {
    let text = await puppeteerPage.$eval(selector, el => el.innerHTML);
    return text;
  },

  async getOneElementInnerText(selector) {
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

  async getAllHistograms(testclass) {
    await puppeteerUtils(puppeteerPage).waitByClass(testclass);
    const histograms = await puppeteerPage.$$eval(
      `[data-testclass=${testclass}]`,
      els => {
        return els.map(el => {
          return el.dataset.testid.substring(
            "histogram_".length,
            el.dataset.testid.length
          );
        });
      }
    );
    return histograms;
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
    const checkboxId = `category-select-${category}`;
    await puppeteerUtils(puppeteerPage).waitByID(checkboxId);
    const checkedPseudoclass = await puppeteerPage.$eval(
      `[data-testid='${checkboxId}']`,
      el => {
        return el.matches(":checked");
      }
    );
    if (!checkedPseudoclass) {
      await puppeteerUtils(puppeteerPage).clickOn(checkboxId);
    }
    try {
      const categoryRow = await puppeteerUtils(puppeteerPage).waitByID(
        `category-expand-${category}`
      );
      const isExpanded = await categoryRow.$(
        "[data-testclass='category-expand-is-expanded']"
      );
      if (isExpanded) {
        await puppeteerUtils(puppeteerPage).clickOn(
          `category-expand-${category}`
        );
      }
    } catch {}
  },

  async calcDragCoordinates(testid, coordinateAsPercent) {
    const el = await puppeteerUtils(puppeteerPage).waitByID(testid);
    const size = await el.boxModel();
    const coords = {
      start: {
        x: Math.floor(size.width * coordinateAsPercent.x1),
        y: Math.floor(size.height * coordinateAsPercent.y1)
      },
      end: {
        x: Math.floor(size.width * coordinateAsPercent.x2),
        y: Math.floor(size.height * coordinateAsPercent.y2)
      }
    };
    return coords;
  },

  async selectCategory(category, values, reset = true) {
    if (reset) await this.resetCategory(category);
    await puppeteerUtils(puppeteerPage).clickOn(`category-expand-${category}`);
    await puppeteerUtils(puppeteerPage).clickOn(`category-select-${category}`);
    for (const val of values) {
      await puppeteerUtils(puppeteerPage).clickOn(
        `categorical-value-select-${category}-${val}`
      );
    }
  },

  async reset() {
    await puppeteerUtils(puppeteerPage).clickOn("reset");
    // loading state never actually happens, reset is too fast
    await page.waitFor(200);
  },

  async errorAppears() {
    await puppeteerUtils(puppeteerPage).waitByClass("toast", {
      timeout: 200
    });
  }
});
