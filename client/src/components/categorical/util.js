// jshint esversion: 6

// values is [ [optVal, optIdx], ...]
// index is range array
// return sorted index

import isNumber from "is-number";
import _ from "lodash";

const sortedCategoryValues = values => {
  /* this sort could be memoized for perf */

  const strings = [];
  const ints = [];

  _.forEach(values, v => {
    if (isNumber(v[0])) {
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

  return ints.concat(strings);
};

export default sortedCategoryValues;
