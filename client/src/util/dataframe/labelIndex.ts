/* eslint-disable max-classes-per-file -- Classes are interrelated*/
/**
Label indexing - map a label to & from an integer offset.  See Dataframe
for how this is used.
**/

import { rangeFill as fillRange } from "../range";
import { __getMemoId } from "./util";

/*
Private utility functions
*/
function extent(tarr: any) {
  let min = 0x7fffffff;
  let max = ~min; // eslint-disable-line no-bitwise -- Establishes 0 of same size
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
  maxOffset: any;

  /*
  identity/noop index, with small assumptions that labels are int32
  */
  constructor(maxOffset: any) {
    this.maxOffset = maxOffset;
  }

  get __id() {
    return `IdentityInt32Index_${this.maxOffset}`;
  }

  labels() {
    // memoize
    const k = fillRange(new Int32Array(this.maxOffset));
    this.labels = function labels() {
      return k;
    };
    return k;
  }

  getOffset(i: any) {
    // label to offset
    return Number.isInteger(i) && i >= 0 && i < this.maxOffset ? i : undefined;
  }

  getOffsets(arr: any) {
    // labels to offsets
    return arr.map((i: any) => this.getOffset(i));
  }

  getLabel(i: any) {
    // offset to label
    return Number.isInteger(i) && i >= 0 && i < this.maxOffset ? i : undefined;
  }

  getLabels(arr: any) {
    // offsets to labels
    return arr.map((i: any) => this.getLabel(i));
  }

  size() {
    return this.maxOffset;
  }

  __promote(labelArray: any) {
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
    // @ts-expect-error ts-migrate(2345) FIXME: Argument of type 'number[]' is not assignable to p... Remove this comment to see the full error message
    return new DenseInt32Index(labelArray, [minLabel, maxLabel]);
  }

  subset(labels: any) {
    /* validate subset */
    const { maxOffset } = this;
    for (let i = 0, l = labels.length; i < l; i += 1) {
      const label = labels[i];
      if (!Number.isInteger(label) || label < 0 || label >= maxOffset)
        throw new RangeError(`offset or label: ${label}`);
    }
    return this.__promote(labels);
  }

  /* identity index - labels are offsets */
  isubset(offsets: any) {
    return this.subset(offsets);
  }

  /* identity index - labels are offsets */
  isubsetMask(mask: any) {
    let count = 0;
    if (mask.length !== this.maxOffset) {
      throw new RangeError("mask has invalid length for index");
    }
    let labels = new Int32Array(mask.length);
    for (let i = 0, l = mask.length; i < l; i += 1) {
      if (mask[i]) {
        labels[count] = i;
        count += 1;
      }
    }
    labels = labels.slice(0, count);
    return this.subset(labels);
  }

  withLabel(label: any) {
    if (label === this.maxOffset) {
      return new IdentityInt32Index(label + 1);
    }
    return this.__promote([...this.labels(), label]);
  }

  withLabels(labels: any) {
    return this.__promote([...this.labels(), ...labels]);
  }

  dropLabel(label: any) {
    if (label === this.maxOffset - 1) {
      return new IdentityInt32Index(label);
    }
    const labelArray = [...this.labels()];
    labelArray.splice(labelArray.indexOf(label), 1);
    return this.__promote(labelArray);
  }
}

class DenseInt32Index {
  __id: any;

  getLabel: any;

  getLabels: any;

  getOffset: any;

  getOffsets: any;

  index: any;

  minLabel: any;

  rindex: any;

  /*
  DenseInt32Index indexes integer labels, and uses Int32Array typed arrays
  for both forward and reverse indexing.   This means that the min/max range
  of the forward index labels must be known a priori (so that the index
  array can be pre-allocated).
  */
  constructor(labels: any, labelRange = null) {
    if (labels.constructor !== Int32Array) {
      labels = new Int32Array(labels);
    }

    if (!labelRange) {
      // @ts-expect-error ts-migrate(2322) FIXME: Type 'number[]' is not assignable to type 'null'.
      labelRange = extent(labels);
    }
    // @ts-expect-error ts-migrate(2488) FIXME: Type 'null' must have a '[Symbol.iterator]()' meth... Remove this comment to see the full error message
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
    this.__id = __getMemoId();
    this.__compile();
  }

  __compile() {
    const { minLabel, index, rindex } = this;
    this.getOffset = function getOffset(l: any) {
      if (!Number.isInteger(l)) return undefined;
      const offset = index[l - minLabel];
      return offset === -1 ? undefined : offset;
    };

    this.getOffsets = function getOffsets(arr: any) {
      return arr.map((i: any) => this.getOffset(i));
    };

    this.getLabel = function getLabel(i: any) {
      return Number.isInteger(i) ? rindex[i] : undefined;
    };

    this.getLabels = function getLabels(arr: any) {
      return arr.map((i: any) => this.getLabel(i));
    };
  }

  labels() {
    return this.rindex;
  }

  size() {
    return this.rindex.length;
  }

  __promote(labelArray: any) {
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
    // @ts-expect-error ts-migrate(2345) FIXME: Argument of type 'number[]' is not assignable to p... Remove this comment to see the full error message
    return new DenseInt32Index(labelArray, [minLabel, maxLabel]);
  }

  subset(labels: any) {
    /* validate subset */
    for (let i = 0, l = labels.length; i < l; i += 1) {
      const label = labels[i];
      const offset = this.getOffset(label);
      if (offset === undefined || offset === -1)
        throw new RangeError(`unknown label: ${label}`);
    }
    return this.__promote(labels);
  }

  isubset(offsets: any) {
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

  isubsetMask(mask: any) {
    const { rindex } = this;
    if (mask.length !== rindex.length)
      throw new RangeError("mask has invalid length for index");
    let count = 0;
    let labels = new Int32Array(mask.length);
    for (let i = 0, l = mask.length; i < l; i += 1) {
      if (mask[i]) {
        labels[count] = rindex[i];
        count += 1;
      }
    }
    labels = labels.slice(0, count);
    return this.__promote(labels);
  }

  withLabel(label: any) {
    return this.__promote([...this.labels(), label]);
  }

  withLabels(labels: any) {
    return this.__promote([...this.labels(), ...labels]);
  }

  dropLabel(label: any) {
    const labelArray = [...this.labels()];
    labelArray.splice(labelArray.indexOf(label), 1);
    return this.__promote(labelArray);
  }
}

class KeyIndex {
  __id: any;

  getLabel: any;

  getLabels: any;

  getOffset: any;

  getOffsets: any;

  index: any;

  rindex: any;

  /*
  KeyIndex indexes arbitrary JS primitive types, and uses a Map()
  as its core data structure.
  */
  constructor(labels: any) {
    const index = new Map();
    if (labels === undefined) {
      labels = [];
    }
    const rindex = labels;
    labels.forEach((v: any, i: any) => {
      index.set(v, i);
    });

    if (index.size !== rindex.length) {
      /* if true, there was a duplicate in the keys */
      throw new Error("duplicate label provided to KeyIndex");
    }

    this.index = index;
    this.rindex = rindex;
    this.__id = __getMemoId();
    this.__compile();
  }

  __compile() {
    const { index, rindex } = this;
    this.getOffset = function getOffset(k: any) {
      return index.get(k);
    };

    this.getOffsets = function getOffsets(arr: any) {
      return arr.map((l: any) => this.getOffset(l));
    };

    this.getLabel = function getLabel(i: any) {
      return Number.isInteger(i) ? rindex[i] : undefined;
    };

    this.getLabels = function getLabels(arr: any) {
      return arr.map((i: any) => this.getLabel(i));
    };
  }

  labels() {
    return this.rindex;
  }

  size() {
    return this.rindex.length;
  }

  subset(labels: any) {
    /* validate subset */
    for (let i = 0, l = labels.length; i < l; i += 1) {
      const label = labels[i];
      const offset = this.getOffset(label);
      if (offset === undefined || offset === -1) {
        throw new RangeError(`unknown label: ${label}`);
      }
    }

    return new KeyIndex(labels);
  }

  isubset(offsets: any) {
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

  isubsetMask(mask: any) {
    const { rindex } = this;
    if (mask.length !== rindex.length)
      throw new RangeError("mask has invalid length for index");
    let labels = new Array(mask.length);
    let count = 0;
    for (let i = 0, l = mask.length; i < l; i += 1) {
      if (mask[i]) {
        labels[count] = rindex[i];
        count += 1;
      }
    }
    labels = labels.slice(0, count);
    return new KeyIndex(labels);
  }

  withLabel(label: any) {
    return new KeyIndex([...this.rindex, label]);
  }

  withLabels(labels: any) {
    return new KeyIndex([...this.rindex, ...labels]);
  }

  dropLabel(label: any) {
    const idx = this.rindex.indexOf(label);
    const labelArray = [...this.rindex];
    labelArray.splice(idx, 1);
    return new KeyIndex(labelArray);
  }
}

function isLabelIndex(i: any) {
  return (
    i instanceof IdentityInt32Index ||
    i instanceof DenseInt32Index ||
    i instanceof KeyIndex
  );
}

export { DenseInt32Index, IdentityInt32Index, KeyIndex, isLabelIndex };
/* eslint-enable max-classes-per-file -- enable*/
