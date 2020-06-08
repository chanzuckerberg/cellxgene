export const puppeteerUtils = (page) => ({
  async waitByID(testId, props = {}) {
    return page.waitForSelector(`[data-testid='${testId}']`, props);
  },

  async waitByClass(testClass, props = {}) {
    return page.waitForSelector(`[data-testclass='${testClass}']`, props);
  },

  async waitForAllByIds(testIds) {
    await Promise.all(
      testIds.map((testId) => page.waitForSelector(`[data-testid='${testId}']`))
    );
  },

  async getAllByClass(testClass) {
    return page.$$eval(`[data-testclass=${testClass}]`, (eles) =>
      eles.map((ele) => ele.dataset.testid)
    );
  },

  async typeInto(testId, text) {
    // blueprint's  typeahead is treating typing weird, clicking & waiting first solves this
    // only works for text without special characters
    await this.waitByID(testId);
    const selector = `[data-testid='${testId}']`;
    // type ahead can be annoying if you don't pause before you type
    await page.click(selector);
    await page.waitFor(200);
    await page.type(selector, text);
  },

  async clearInputAndTypeInto(testId, text) {
    await this.waitByID(testId);
    const selector = `[data-testid='${testId}']`;
    // only works for text without special characters
    // type ahead can be annoying if you don't pause before you type
    await page.click(selector);
    await page.waitFor(200);
    // select all
    await page.click(selector, { clickCount: 3 });
    await page.keyboard.press("Backspace");
    await page.type(selector, text);
  },

  async clickOn(testid, options = {}) {
    await this.waitByID(testid, options);
    const click = await page.click(`[data-testid='${testid}']`);
    await page.waitFor(50);
    return click;
  },

  async getOneElementInnerHTML(selector, options = {}) {
    await page.waitForSelector(selector, options);
    return page.$eval(selector, (el) => el.innerHTML);
  },

  async getOneElementInnerText(selector) {
    await page.waitForSelector(selector);
    return page.$eval(selector, (el) => el.innerText);
  },

  async getElementCoordinates(testid) {
    return page.$eval(`[data-testid='${testid}']`, (elem) => {
      const { left, top } = elem.getBoundingClientRect();
      return [left, top];
    });
  },
});
