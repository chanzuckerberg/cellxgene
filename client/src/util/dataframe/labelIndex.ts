/* eslint-disable max-classes-per-file -- Classes are interrelated*/
/**
Label indexing - map a label to & from an integer offset.  See Dataframe
for how this is used.
**/

import { AnyArray, TypedArray } from "../../common/types/arraytypes";
import { rangeFill as fillRange } from "../range";
import { __getMemoId } from "./util";

export type OffsetArray =
  | Int8Array
  | Uint8Array
  | Int16Array
  | Uint16Array
  | Int32Array
  | Uint32Array
  | number[];

/*
Private utility functions
*/

/** @internal */
function extent(tarr: OffsetArray): [number, number] {
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

abstract class LabelIndex<LabelArrayType extends AnyArray> {
  readonly __id: string;

  constructor(id: string) {
    this.__id = id;
  }

  abstract labels(): LabelArrayType;

  abstract getOffset(label: LabelArrayType[0]): number | undefined;

  abstract getOffsets(labels: LabelArrayType): (number | undefined)[];

  abstract getLabel(offset: number): LabelArrayType[0] | undefined;

  abstract getLabels(
    offsets: OffsetArray
  ): Array<LabelArrayType[0] | undefined>;

  abstract size(): number;

  abstract subset(labels: LabelArrayType): LabelIndexImplementations;

  abstract isubset(offsets: OffsetArray): LabelIndexImplementations;

  abstract isubsetMask(mask: Uint8Array): LabelIndexImplementations;

  abstract withLabel(label: LabelArrayType[0]): LabelIndexImplementations;

  abstract withLabels(labels: LabelArrayType): LabelIndexImplementations;

  abstract dropLabel(label: LabelArrayType[0]): LabelIndexImplementations;
}

class IdentityInt32Index extends LabelIndex<Int32Array> {
  readonly maxOffset: Int32Array[0];

  /*
  identity/noop index, with small assumptions that labels are int32
  */
  constructor(maxOffset: number) {
    super(`IdentityInt32Index_${maxOffset}`);
    this.maxOffset = maxOffset;
  }

  labels(): Int32Array {
    // memoize
    const k = fillRange(new Int32Array(this.maxOffset));
    this.labels = function labels() {
      return k;
    };
    return k;
  }

  getOffset(label: number): number | undefined {
    // label to offset
    return Number.isInteger(label) && label >= 0 && label < this.maxOffset
      ? label
      : undefined;
  }

  getOffsets(labels: OffsetArray): (number | undefined)[] {
    // labels to offsets
    const result = new Array(labels.length);
    for (let i = 0; i < labels.length; i += 1) {
      result[i] = this.getOffset(labels[i]);
    }
    return result;
  }

  getLabel(offset: number): number | undefined {
    // offset to label
    return Number.isInteger(offset) && offset >= 0 && offset < this.maxOffset
      ? offset
      : undefined;
  }

  getLabels(offsets: OffsetArray): (number | undefined)[] {
    // offsets to labels
    const result = new Array(offsets.length);
    for (let i = 0; i < offsets.length; i += 1) {
      result[i] = this.getOffset(offsets[i]);
    }
    return result;
  }

  size(): number {
    return this.maxOffset;
  }

  __promote(labelArray: OffsetArray): LabelIndexImplementations {
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

  subset(labels: OffsetArray): LabelIndexImplementations {
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
  isubset(offsets: OffsetArray): LabelIndexImplementations {
    return this.subset(offsets);
  }

  /* identity index - labels are offsets */
  isubsetMask(mask: Uint8Array): LabelIndexImplementations {
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

  withLabel(label: number): LabelIndexImplementations {
    if (label === this.maxOffset) {
      return new IdentityInt32Index(label + 1);
    }
    return this.__promote([...this.labels(), label]);
  }

  withLabels(labels: OffsetArray): LabelIndexImplementations {
    return this.__promote([...this.labels(), ...labels]);
  }

  dropLabel(label: number): LabelIndexImplementations {
    if (label === this.maxOffset - 1) {
      return new IdentityInt32Index(label);
    }
    const labelArray = [...this.labels()];
    labelArray.splice(labelArray.indexOf(label), 1);
    return this.__promote(labelArray);
  }
}

class DenseInt32Index extends LabelIndex<Int32Array> {
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
    super(__getMemoId());
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

class KeyIndex extends LabelIndex {
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
    super();
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

  size(): number {
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

type LabelIndexImplementations =
  | DenseInt32Index
  | IdentityInt32Index
  | KeyIndex;

function isLabelIndex(i: unknown): i is LabelIndex {
  return (
    i instanceof IdentityInt32Index ||
    i instanceof DenseInt32Index ||
    i instanceof KeyIndex
  );
}

export { DenseInt32Index, IdentityInt32Index, KeyIndex, isLabelIndex };
/* eslint-enable max-classes-per-file -- enable*/
