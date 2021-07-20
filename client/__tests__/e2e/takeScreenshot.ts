function toFilename(name: any) {
  return name.replace(/[^a-z0-9.-]+/gi, "-");
}

// @ts-expect-error ts-migrate(2451) FIXME: Cannot redeclare block-scoped variable 'takeScreen... Remove this comment to see the full error message
async function takeScreenshot(currentTestName: any, page: any) {
  const testName = toFilename(currentTestName);

  // Take a screenshot at the point of failure
  const date = new Date().toISOString();
  const screenshotName = `${date}-${testName}.png`;

  await page.screenshot({
    path: `./__tests__/screenshots/ignoreE2E-screenshot-${screenshotName}`,
  });
}

module.exports = takeScreenshot;
