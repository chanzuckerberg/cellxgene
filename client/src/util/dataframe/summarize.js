/*
Private dataframe support functions
*/

export function summarizeContinuous(col) {
  let min;
  let max;
  let nan = 0;
  let pinf = 0;
  let ninf = 0;
  if (col) {
    for (let r = 0, l = col.length; r < l; r += 1) {
      const val = Number(col[r]);
      if (Number.isFinite(val)) {
        if (min === undefined) {
          min = val;
          max = val;
        } else {
          min = val < min ? val : min;
          max = val > max ? val : max;
        }
      } else if (Number.isNaN(val)) {
        nan += 1;
      } else if (val > 0) {
        pinf += 1;
      } else {
        ninf += 1;
      }
    }
  }
  return {
    categorical: false,
    min,
    max,
    nan,
    pinf,
    ninf
  };
}

export function summarizeCategorical(col) {
  const categoryCounts = new Map();
  if (col) {
    for (let r = 0, l = col.length; r < l; r += 1) {
      const val = col[r];
      let curCount = categoryCounts.get(val);
      if (curCount === undefined) curCount = 0;
      categoryCounts.set(val, curCount + 1);
    }
  }
  return {
    categorical: true,
    categories: [...categoryCounts.keys()],
    categoryCounts,
    numCategories: categoryCounts.size
  };
}
