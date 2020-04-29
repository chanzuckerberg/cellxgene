const cheerio = require("cheerio");
const crypto = require("crypto");
HtmlWebpackPlugin = require("html-webpack-plugin");

class CspHashPlugin {
  constructor(opts) {
    this.opts = { ...opts };
  }

  apply(compiler) {
    compiler.hooks.compilation.tap("CspHashPlugin", (compilation) => {
      HtmlWebpackPlugin.getHooks(compilation).beforeEmit.tapAsync(
        "CspHashPlugin",
        (data, cb) => {
          const { filename } = this.opts;

          if (filename) {
            const $ = cheerio.load(data.html, { decodeEntities: false });
            const results = {};
            results["script-hashes"] = $("script:not([src])")
              .map((i, elmt) => this.digest($(elmt).html()))
              .get();
            results["style-hashes"] = $("style:not([href])")
              .map((i, elmt) => this.digest($(elmt).html()))
              .get();

            const json = JSON.stringify(results);
            compilation.assets[filename] = {
              source: () => json,
              size: () => json.length,
            };
          }

          // Tell webpack to move on
          cb(null, data);
        }
      );
    });
  }

  digest(str) {
    const hash = crypto
      .createHash("sha256")
      .update(str, "utf8")
      .digest("base64");
    return `sha256-${  hash}`;
  }
}

module.exports = CspHashPlugin;
