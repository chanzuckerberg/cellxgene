export const puppeteerUtils = puppeteerPage => ({
  async waitByID(testid) {
    return await puppeteerPage.waitForSelector(`[data-testid='${testid}']`);
  },

  async waitByClass(testclass) {
    return await puppeteerPage.waitForSelector(
      `[data-testclass='${testclass}']`
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
  async drag(el_box, start, end, lasso = false) {
    const x1 = el_box.content[0].x + start.x;
    const x2 = el_box.content[0].x + end.x;
    const y1 = el_box.content[0].y + start.y;
    const y2 = el_box.content[0].y + end.y;
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
      divs => {
        return divs.map(div => {
          // TODO get from ids
          return div.id.substring("histogram_".length, div.id.length);
        });
      }
    );
    return histograms;
  },

  // make sure these are children of the category
  async getAllCategories(category) {
    await c.waitByClass("categorical-row");
    const categories = await puppeteerPage.$$eval(
      `[data-testid="category-${category}"] [data-testclass=categorical-row]`,
      divs => {
        return divs.map(div => {
          // TODO get from ids
          const val = div.querySelector("[data-testclass='categorical-value']");
          return val.dataset.testid.substring(
            "categorical-value-".length,
            val.dataset.testid.length
          );
        });
      }
    );
    return categories;
  },

  async getAllCategoriesAndCounts(category) {
    await puppeteerUtils(puppeteerPage).waitByClass("categorical-row");
    const categories = await puppeteerPage.$$eval(
      `[data-testid="category-${category}"] [data-testclass=categorical-row]`,
      divs => {
        let result = {};
        divs.forEach(div => {
          const cat = div.querySelector("[data-testclass='categorical-value']")
            .innerText;
          const count = div.querySelector(
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
    await puppeteerUtils(puppeteerPage).clickOn(
      `cellset-button-${num}`
    );
    return await puppeteerUtils(puppeteerPage).getOneElementInnerText(
       `[data-testid='cellset-count-${num}']`
    );
  },

  async reset() {
    await puppeteerUtils(puppeteerPage).clickOn("reset");
    await page.waitFor(100)
  }
});

// TODO cellxgene specific
// find histograms of type
// reset page
// select cells

// TODO fixtures for different datasets
