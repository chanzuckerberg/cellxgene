/*
    via https://github.com/nodef/extra-number/blob/master/scripts/significantDigits.js
*/

const significantDigits = (n: number): number =>
  n
    .toExponential()
    .replace(/e[+\-0-9]*$/, "")
    .replace(/^0\.?0*|\./, "").length;

export default significantDigits;
