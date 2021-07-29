/* eslint-disable import/no-extraneous-dependencies -- this file is a devDependency*/
/* eslint-disable @blueprintjs/classes-constants -- we don't import blueprint here  */
// eslint-disable-next-line @typescript-eslint/no-var-requires --- FIXME: disabled temporarily on migrate to TS.
const cheerio = require("cheerio");
// eslint-disable-next-line @typescript-eslint/no-var-requires --- FIXME: disabled temporarily on migrate to TS.
const crypto = require("crypto");
// eslint-disable-next-line @typescript-eslint/no-var-requires --- FIXME: disabled temporarily on migrate to TS.
const HtmlWebpackPlugin = require("html-webpack-plugin");

const digest = (str) => {
  const hash = crypto.createHash("sha256").update(str, "utf8").digest("base64");
  return `sha256-${hash}`;
};
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

          const $ = cheerio.load(data.html, { decodeEntities: false });

          if (filename) {
            const results = {};
            results["script-hashes"] = $("script:not([src]):not([no-csp-hash])")
              .map((i, elmt) => digest($(elmt).html()))
              .get();

            const json = JSON.stringify(results);
            compilation.assets[filename] = {
              source: () => json,
              size: () => json.length,
            };
          }

          // Remove no-csp-hash attributes. Cheerio does not parse Jinja templates
          // correctly, so we brute force this with a regular expression.
          data.html = data.html.replace(
            /(<script .*)no-csp-hash(.*>)/,
            (match, p1, p2) => [p1, p2].join("")
          );

          // Tell webpack to move on
          cb(null, data);
        }
      );
    });
  }
}

module.exports = CspHashPlugin;
/* eslint-enable @blueprintjs/classes-constants -- we don't import blueprint here  */
/* eslint-enable import/no-extraneous-dependencies -- enable*/
