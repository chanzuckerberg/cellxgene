/* 
    via https://github.com/nodef/extra-number/blob/master/scripts/significantDigits.js
*/

// eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any -- - FIXME: disabled temporarily on migrate to TS.
export default (n: any) =>
  n
    .toExponential()
    .replace(/e[+\-0-9]*$/, "")
    .replace(/^0\.?0*|\./, "").length;
