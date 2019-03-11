export const puppeteerUtils = puppeteerPage => ({
  async getOneElementInnerHTML(selector) {
    let text = await puppeteerPage.$eval(selector, el => el.innerHTML);
    return text;
  },

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
  }
});
