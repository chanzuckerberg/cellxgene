const { setup } = require("jest-environment-puppeteer");

module.exports = async () => {
  console.log("Global setup...");
  await new Promise(setTimeout(setup, 30000));
};
