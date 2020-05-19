/* 
    via https://github.com/nodef/extra-number/blob/master/scripts/significantDigits.js
*/

export default (n) => {
  return n
    .toExponential()
    .replace(/e[+\-0-9]*$/, "")
    .replace(/^0\.?0*|\./, "").length;
};
