/*
Sort category values (labels) in the order we want for presentation.

TL;DR: sort order is:
* numbers or number-like strings first, in numeric order
* most strings, in case-insenstive unicode sort order
* then 'nan' (any case)
* then, IF isUseAnno is true, globals.unassignedCategoryLabel
*/

import isNumber from "is-number";
import * as globals from "../globals";

function caseInsensitiveCompare(a: any, b: any) {
  const textA = String(a).toUpperCase();
  const textB = String(b).toUpperCase();
  return textA < textB ? -1 : textA > textB ? 1 : 0;
}

const catLabelSort = (isUserAnno: any, values: any) => {
  /* this sort could be memoized for perf */

  const strings: any = [];
  const ints: any = [];
  const unassignedOrNaN: any = [];

  values.forEach((v: any) => {
    if (isUserAnno && v === globals.unassignedCategoryLabel) {
      unassignedOrNaN.push(v);
    } else if (String(v).toLowerCase() === "nan") {
      unassignedOrNaN.push(v);
    } else if (isNumber(v)) {
      ints.push(v);
    } else {
      strings.push(v);
    }
  });

  strings.sort(caseInsensitiveCompare);
  // @ts-expect-error ts-migrate(7006) FIXME: Parameter 'a' implicitly has an 'any' type.
  ints.sort((a, b) => +a - +b);
  unassignedOrNaN.sort(caseInsensitiveCompare);

  return ints.concat(strings, unassignedOrNaN);
};

export default catLabelSort;
