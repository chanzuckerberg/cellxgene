/* eslint-disable max-classes-per-file */
/**
Label indexing - map a label to & from an integer offset.  See Dataframe
for how this is used.
**/

import { rangeFill as fillRange } from "../range";

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

class IdentityInt32Index {
  /*
  identity/noop index, with small assumptions that labels are int32
  */
  constructor(maxOffset) {
    this.maxOffset = maxOffset;
  }

  labels() {
    // memoize
    const k = fillRange(new Int32Array(this.maxOffset));
    this.labels = function labels() {
      return k;
    };
    return k;
  }

  // eslint-disable-next-line class-methods-use-this
  getOffset(i) {
    // label to offset
    return i;
  }

  // eslint-disable-next-line class-methods-use-this
  getOffsets(arr) {
    // labels to offsets
    return arr;
  }

  // eslint-disable-next-line class-methods-use-this
  getLabel(i) {
    // offset to label
    return i;
  }

  // eslint-disable-next-line class-methods-use-this
  getLabels(arr) {
    // offsets to labels
    return arr;
  }

  size() {
    return this.maxOffset;
  }

  __promote(labelArray) {
    /*
    time/space decision - based on the resulting density
    */
    const [minLabel, maxLabel] = extent(labelArray);
    if (minLabel === 0 && maxLabel === labelArray.length - 1)
      return new IdentityInt32Index(labelArray.length);

    const labelSpaceSize = maxLabel - minLabel + 1;
    const density = labelSpaceSize / this.maxOffset;
    /* 0.1 is a magic number, that needs testing to optimize */
    if (density < 0.1) {
      return new KeyIndex(labelArray);
    }
    return new DenseInt32Index(labelArray, [minLabel, maxLabel]);
  }

  subset(labels) {
    /* validate subset */
    const { maxOffset } = this;
    for (let i = 0, l = labels.length; i < l; i += 1) {
      const label = labels[i];
      if (label < 0 || label >= maxOffset)
        throw new RangeError(`offset or label: ${label}`);
    }
    return this.__promote(labels);
  }

  /* identity index - labels are offsets */
  isubset(offsets) {
    return this.subset(offsets);
  }

  withLabel(label) {
    if (label === this.maxOffset) {
      return new IdentityInt32Index(label + 1);
    }
    return this.__promote([...this.labels(), label]);
  }

  withLabels(labels) {
    return this.__promote([...this.labels(), ...labels]);
  }

  dropLabel(label) {
    if (label === this.maxOffset - 1) {
      return new IdentityInt32Index(label);
    }
    const labelArray = [...this.labels()];
    labelArray.splice(labelArray.indexOf(label), 1);
    return this.__promote(labelArray);
  }
}

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

    this.getOffsets = function getOffsets(arr) {
      const res = new arr.constructor(arr.length);
      for (let i = 0, len = arr.length; i < len; i += 1) {
        res[i] = index[arr[i] - minLabel];
      }
      return res;
    };

    this.getLabel = function getLabel(i) {
      return rindex[i];
    };

    this.getLabels = function getLabels(arr) {
      const res = new arr.constructor(arr.length);
      for (let i = 0, len = arr.length; i < len; i += 1) {
        res[i] = rindex[arr[i]];
      }
      return res;
    };
  }

  labels() {
    return this.rindex;
  }

  size() {
    return this.rindex.length;
  }

  __promote(labelArray) {
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

  subset(labels) {
    /* validate subset */
    for (let i = 0, l = labels.length; i < l; i += 1) {
      const label = labels[i];
      const offset = this.getOffset(label);
      if (offset === undefined || offset === -1)
        throw new RangeError(`unknown label: ${label}`);
    }

    return this.__promote(labels);
  }

  // eslint-disable-next-line class-methods-use-this
  isubset(offsets) {
    /* validate subset */
    const { rindex } = this;
    const maxOffset = rindex.length;
    const labels = new Int32Array(offsets.length);
    for (let i = 0, l = offsets.length; i < l; i += 1) {
      const offset = offsets[i];
      if (offset < 0 || offset >= maxOffset)
        throw new RangeError(`out of bounds offset: ${offset}`);
      labels[i] = rindex[offset];
    }

    return this.__promote(labels);
  }

  withLabel(label) {
    return this.__promote([...this.labels(), label]);
  }

  withLabels(labels) {
    return this.__promote([...this.labels(), ...labels]);
  }

  dropLabel(label) {
    const labelArray = [...this.labels()];
    labelArray.splice(labelArray.indexOf(label), 1);
    return this.__promote(labelArray);
  }
}

class KeyIndex {
  /*
  KeyIndex indexes arbitrary JS primitive types, and uses a Map()
  as its core data structure.
  */
  constructor(labels) {
    const index = new Map();
    if (labels === undefined) {
      labels = [];
    }
    const rindex = labels;
    labels.forEach((v, i) => {
      index.set(v, i);
    });

    if (index.size !== rindex.length) {
      /* if true, there was a duplicate in the keys */
      throw new Error("duplicate label provided to KeyIndex");
    }

    this.index = index;
    this.rindex = rindex;
    this.__compile();
  }

  __compile() {
    const { index, rindex } = this;
    this.getOffset = function getOffset(k) {
      return index.get(k);
    };

    this.getOffsets = function getOffsets(arr) {
      const res = new arr.constructor(arr.length);
      for (let i = 0, len = arr.length; i < len; i += 1) {
        res[i] = index.get(arr[i]);
      }
      return res;
    };

    this.getLabel = function getLabel(i) {
      return rindex[i];
    };

    this.getLabels = function getLabels(arr) {
      const res = new arr.constructor(arr.length);
      for (let i = 0, len = arr.length; i < len; i += 1) {
        res[i] = rindex[arr[i]];
      }
      return res;
    };
  }

  labels() {
    return this.rindex;
  }

  size() {
    return this.rindex.length;
  }

  subset(labels) {
    /* validate subset */
    for (let i = 0, l = labels.length; i < l; i += 1) {
      const label = labels[i];
      const offset = this.getOffset(label);
      if (offset === undefined) throw new RangeError(`unknown label: ${label}`);
    }

    return new KeyIndex(labels);
  }

  // eslint-disable-next-line class-methods-use-this
  isubset(offsets) {
    const { rindex } = this;
    const maxOffset = rindex.length;
    const labels = new Array(offsets.length);
    for (let i = 0, l = offsets.length; i < l; i += 1) {
      const offset = offsets[i];
      if (offset < 0 || offset >= maxOffset)
        throw new RangeError(`out of bounds offset: ${offset}`);
      labels[i] = rindex[offset];
    }

    return new KeyIndex(labels);
  }

  withLabel(label) {
    return new KeyIndex([...this.rindex, label]);
  }

  withLabels(labels) {
    return new KeyIndex([...this.rindex, ...labels]);
  }

  dropLabel(label) {
    const idx = this.rindex.indexOf(label);
    const labelArray = [...this.rindex];
    labelArray.splice(idx, 1);
    return new KeyIndex(labelArray);
  }
}
/* eslint-enable class-methods-use-this */

function isLabelIndex(i) {
  return (
    i instanceof IdentityInt32Index ||
    i instanceof DenseInt32Index ||
    i instanceof KeyIndex
  );
}

export { DenseInt32Index, IdentityInt32Index, KeyIndex, isLabelIndex };
