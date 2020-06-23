function toFilename(name) {
  return name.replace(/[^a-z0-9.-]+/gi, "-");
}

async function takeScreenshot(currentTestName, page) {
  const testName = toFilename(currentTestName);

  // Take a screenshot at the point of failure
  const date = new Date().toISOString();
  const screenshotName = `${date}-${testName}.png`;

  await page.screenshot({
    path: `./__tests__/screenshots/ignoreE2E-screenshot-${screenshotName}`,
  });
}

module.exports = takeScreenshot;
