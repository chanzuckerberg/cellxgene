// jshint esversion: 6

// values is [ [optVal, optIdx], ...]
// index is range array
// return sorted index

/*
Sort category values (labels) in the order we want for presentation.

TL;DR: numeric sort for number-like strings, then strings in case-ignoring alpha
order.  Except, when isUserAnno is true, pin globals.unassignedCategoryLabel 
to the end.

*/

import isNumber from "is-number";
import * as globals from "../../globals";

const sortedCategoryValues = (isUserAnno, values) => {
  /* this sort could be memoized for perf */

  const strings = [];
  const ints = [];
  const unassigned = [];

  values.forEach(v => {
    if (isUserAnno && v[0] === globals.unassignedCategoryLabel) {
      unassigned.push(v);
    } else if (isNumber(v[0])) {
      ints.push(v);
    } else {
      strings.push(v);
    }
  });

  strings.sort((a, b) => {
    const textA = String(a[0]).toUpperCase();
    const textB = String(b[0]).toUpperCase();
    return textA < textB ? -1 : textA > textB ? 1 : 0;
  });

  ints.sort((a, b) => +a[0] - +b[0]);

  return ints.concat(strings, unassigned);
};

export default sortedCategoryValues;
