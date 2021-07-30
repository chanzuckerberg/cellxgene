/* eslint-disable no-await-in-loop -- await in loop is needed to emulate sequential user actions */
// eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
export function getTestId(id: any) {
  return `[data-testid='${id}']`;
}

// eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
export function getTestClass(className: any) {
  return `[data-testclass='${className}']`;
}

// eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
export async function waitByID(testId: any, props = {}) {
  return page.waitForSelector(getTestId(testId), props);
}

// eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
export async function waitByClass(testClass: any, props = {}) {
  return page.waitForSelector(`[data-testclass='${testClass}']`, props);
}

// eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
export async function waitForAllByIds(testIds: any) {
  await Promise.all(
    // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
    testIds.map((testId: any) => page.waitForSelector(getTestId(testId)))
  );
}

// eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
export async function getAllByClass(testClass: any) {
  return page.$$(`[data-testclass=${testClass}]`);
}

// eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
export async function typeInto(testId: any, text: any) {
  // blueprint's  typeahead is treating typing weird, clicking & waiting first solves this
  // only works for text without special characters
  await waitByID(testId);
  const selector = getTestId(testId);
  // type ahead can be annoying if you don't pause before you type
  await page.click(selector);
  await page.waitForTimeout(200);
  await page.type(selector, text);
}

// eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
export async function clearInputAndTypeInto(testId: any, text: any) {
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

// eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
export async function clickOn(testId: any, options = {}) {
  await expect(page).toClick(getTestId(testId), options);
}

/**
 * (thuang): There are times when Puppeteer clicks on a button and the page doesn't respond.
 * So I added clickOnUntil() to retry clicking until a given condition is met.
 */
// eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
export async function clickOnUntil(testId: any, assert: any) {
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

// eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
export async function getOneElementInnerHTML(selector: any, options = {}) {
  await page.waitForSelector(selector, options);

  return page.$eval(selector, (el) => el.innerHTML);
}

// eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
export async function getOneElementInnerText(selector: any) {
  expect(page).toMatchElement(selector);

  // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
  return page.$eval(selector, (el) => (el as any).innerText);
}

// eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
export async function getElementCoordinates(testId: any) {
  return page.$eval(getTestId(testId), (elem) => {
    const { left, top } = elem.getBoundingClientRect();
    return [left, top];
  });
}

async function clickTermsOfService() {
  // @ts-expect-error ts-migrate(2554) FIXME: Expected 2 arguments, but got 1.
  if (!(await isElementPresent(getTestId("tos-cookies-accept")))) return;

  await clickOn("tos-cookies-accept");
}

async function nameNewAnnotation() {
  // @ts-expect-error ts-migrate(2554) FIXME: Expected 2 arguments, but got 1.
  if (await isElementPresent(getTestId("annotation-dialog"))) {
    await typeInto("new-annotation-name", "ignoreE2E");
    await clickOn("submit-annotation");

    // wait for the page to load
    await waitByClass("autosave-complete");
  }
}

// eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
export async function goToPage(url: any) {
  await page.goto(url, {
    waitUntil: "networkidle0",
  });

  await nameNewAnnotation();
  await clickTermsOfService();
}

// eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
export async function isElementPresent(selector: any, options: any) {
  // @ts-expect-error ts-migrate(2554) FIXME: Expected 1 arguments, but got 2.
  return Boolean(await page.$(selector, options));
}
/* eslint-enable no-await-in-loop -- await in loop is needed to emulate sequential user actions */
