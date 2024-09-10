const { setup } = require("jest-environment-puppeteer");

module.exports = async () => {
  console.log("Global setup...");
  jest.setTimeout(20 * 1000);
  await setup();
};
