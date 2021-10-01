/* eslint-disable no-await-in-loop -- await in loop is needed to emulate sequential user actions */
export function getTestId(id) {
  return `[data-testid='${id}']`;
}

export function getTestClass(className) {
  return `[data-testclass='${className}']`;
}

export async function waitByID(testId, props = {}) {
  return page.waitForSelector(getTestId(testId), props);
}

export async function waitByClass(testClass, props = {}) {
  return page.waitForSelector(`[data-testclass='${testClass}']`, props);
}

export async function waitForAllByIds(testIds) {
  await Promise.all(
    testIds.map((testId) => page.waitForSelector(getTestId(testId)))
  );
}

export async function getAllByClass(testClass) {
  return page.$$(`[data-testclass=${testClass}]`);
}

export async function typeInto(testId, text) {
  // blueprint's  typeahead is treating typing weird, clicking & waiting first solves this
  // only works for text without special characters
  await waitByID(testId);
  const selector = getTestId(testId);
  // type ahead can be annoying if you don't pause before you type
  await page.click(selector);
  await page.waitForTimeout(200);
  await page.type(selector, text);
}

export async function clearInputAndTypeInto(testId, text) {
  await waitByID(testId);
  const selector = getTestId(testId);
  // only works for text without special characters
  // type ahead can be annoying if you don't pause before you type
  await page.click(selector);
  await page.waitForTimeout(200);
  // select all
  await page.click(selector, { clickCount: 3 });
  await page.keyboard.press("Backspace");
  await page.type(selector, text);
}

export async function clickOn(testId, options = {}) {
  await expect(page).toClick(getTestId(testId), options);
}

/**
 * (thuang): There are times when Puppeteer clicks on a button and the page doesn't respond.
 * So I added clickOnUntil() to retry clicking until a given condition is met.
 */
export async function clickOnUntil(testId, assert) {
  const MAX_RETRY = 10;
  const WAIT_FOR_MS = 200;

  let retry = 0;

  while (retry < MAX_RETRY) {
    try {
      await clickOn(testId);
      await assert();

      break;
    } catch (error) {
      retry += 1;

      await page.waitForTimeout(WAIT_FOR_MS);
    }
  }

  if (retry === MAX_RETRY) {
    throw Error("clickOnUntil() assertion failed!");
  }
}

export async function getOneElementInnerHTML(selector, options = {}) {
  await page.waitForSelector(selector, options);

  return page.$eval(selector, (el) => el.innerHTML);
}

export async function getOneElementInnerText(selector) {
  expect(page).toMatchElement(selector);

  return page.$eval(selector, (el) => el.innerText);
}

export async function getElementCoordinates(testId) {
  return page.$eval(getTestId(testId), (elem) => {
    const { left, top } = elem.getBoundingClientRect();
    return [left, top];
  });
}

async function nameNewAnnotation() {
  if (await isElementPresent(getTestId("annotation-dialog"))) {
    await typeInto("new-annotation-name", "ignoreE2E");
    await clickOn("submit-annotation");

    // wait for the page to load
    await waitByClass("autosave-complete");
  }
}

export async function goToPage(url) {
  await page.goto(url, {
    waitUntil: "networkidle0",
  });

  await nameNewAnnotation();
}

export async function isElementPresent(selector, options) {
  return Boolean(await page.$(selector, options));
}
/* eslint-enable no-await-in-loop -- await in loop is needed to emulate sequential user actions */
