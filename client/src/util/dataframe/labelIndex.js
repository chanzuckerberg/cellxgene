/**
Label indexing - map a label to & from an integer offset.  See Dataframe
for how this is used.
**/

/*
Private utility functions
*/
function extent(tarr) {
  let min = 0x7fffffff;
  let max = ~min; // eslint-disable-line no-bitwise
  for (let i = 0, l = tarr.length; i < l; i += 1) {
    const v = tarr[i];
    if (v < min) {
      min = v;
    }
    if (v > max) {
      max = v;
    }
  }
  return [min, max];
}

function fillRange(arr, start = 0) {
  const larr = arr;
  for (let i = 0, l = larr.length; i < l; i += 1) {
    larr[i] = i + start;
  }
  return larr;
}

/* eslint-disable class-methods-use-this */
class IdentityInt32Index {
  /*
  identity/noop index, with small assumptions that labels are int32
  */
  constructor(maxOffset) {
    this.maxOffset = maxOffset;
  }

  keys() {
    // memoize
    const k = fillRange(new Int32Array(this.maxOffset));
    this.keys = function keys() {
      return k;
    };
    return k;
  }

  getOffset(i) {
    // label to offset
    return i;
  }

  getLabel(i) {
    // offset to label
    return i;
  }

  getMaxOffset() {
    return this.maxOffset;
  }

  cut(labelArray) {
    /*
    if density of resulting integer
    */
    const [minLabel, maxLabel] = extent(labelArray);
    const labelSpaceSize = maxLabel - minLabel + 1;
    const density = labelSpaceSize / this.maxOffset;
    /* 0.1 is a magic number, that needs testing to optimize */
    if (density < 0.1) {
      return new KeyIndex(labelArray);
    }
    return new DenseInt32Index(labelArray, [minLabel, maxLabel]);
  }
}
/* eslint-enable class-methods-use-this */

/* eslint-disable class-methods-use-this */
class DenseInt32Index {
  /*
  DenseInt32Index indexes integer labels, and uses Int32Array typed arrays
  for both forward and reverse indexing.   This means that the min/max range
  of the forward index labels must be known a priori (so that the index
  array can be pre-allocated).
  */
  constructor(labels, labelRange = null) {
    if (labels.constructor !== Int32Array) {
      labels = new Int32Array(labels);
    }

    if (!labelRange) {
      labelRange = extent(labels);
    }
    const [minLabel, maxLabel] = labelRange;
    const labelSpaceSize = maxLabel - minLabel + 1;
    const index = new Int32Array(labelSpaceSize).fill(-1);
    for (let i = 0, l = labels.length; i < l; i += 1) {
      const label = labels[i];
      index[label - minLabel] = i;
    }

    this.minLabel = minLabel;
    this.rindex = labels;
    this.index = index;
    this.__compile();
  }

  __compile() {
    const { minLabel, index, rindex } = this;
    this.getOffset = function getOffset(l) {
      return index[l - minLabel];
    };
    this.getLabel = function getLabel(i) {
      return rindex[i];
    };
  }

  keys() {
    return this.rindex;
  }

  getMaxOffset() {
    return this.rindex.length;
  }

  cut(labelArray) {
    /*
    time/space decision - if we are going to use less than 10% of the
    dense index space, switch to a KeyIndex (which is slower, but uses
    less memory for sparse label spaces).
    */
    const [minLabel, maxLabel] = extent(labelArray);
    const labelSpaceSize = maxLabel - minLabel + 1;
    const density = labelSpaceSize / this.rindex.length;
    /* 0.1 is a magic number, that needs testing to optimize */
    if (density < 0.1) {
      return new KeyIndex(labelArray);
    }
    return new DenseInt32Index(labelArray, [minLabel, maxLabel]);
  }
}
/* eslint-enable class-methods-use-this */

/* eslint-disable class-methods-use-this */
class KeyIndex {
  /*
  KeyIndex indexes arbitrary JS primitive types, and uses a Map()
  as its core data structure.
  */
  constructor(labels) {
    const index = new Map();
    const rindex = labels;
    labels.forEach((v, i) => {
      index.set(v, i);
    });

    this.index = index;
    this.rindex = rindex;
    this.__compile();
  }

  __compile() {
    const { index, rindex } = this;
    this.getOffset = function getOffset(k) {
      return index.get(k);
    };
    this.getLabel = function getLabel(i) {
      return rindex[i];
    };
  }

  keys() {
    return this.rindex;
  }

  getMaxOffset() {
    return this.rindex.length;
  }

  cut(labelArray) {
    return new KeyIndex(labelArray);
  }
}
/* eslint-enable class-methods-use-this */

export { DenseInt32Index, IdentityInt32Index, KeyIndex };
