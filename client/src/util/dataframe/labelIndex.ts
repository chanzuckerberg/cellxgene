/* eslint-disable max-classes-per-file -- Classes are interrelated*/

/**
Label indexing - map a label to & from an integer offset.  See Dataframe
for how this is used.
**/

import { rangeFill as fillRange } from "../range";
import { __getMemoId } from "./util";
import { OffsetArray, LabelType, LabelArray, GenericLabelArray } from "./types";

export abstract class LabelIndexBase {
  readonly __id: string; // memoization helper

  constructor(id: string) {
    this.__id = id;
  }

  abstract labels(): LabelArray;

  /**
   * Look up the offset for the label.
   *
   * @param label - label to look up
   * @returns - offset number or -1 if not found.
   */
  abstract getOffset(label: LabelType): number;

  getOffsets(labels: LabelArray): Int32Array {
    // labels to offsets
    const result = new Int32Array(labels.length);
    for (let i = 0; i < labels.length; i += 1) {
      result[i] = this.getOffset(labels[i]);
    }
    return result;
  }

  /**
   * Look up the label for the offset.
   *
   * @param offset - offset to look up
   * @returns - label or undefined if not found.
   */
  abstract getLabel(offset: number): LabelType | undefined;

  getLabels(offsets: OffsetArray): (LabelType | undefined)[] {
    // offsets to labels
    const result = new Array(offsets.length);
    for (let i = 0; i < offsets.length; i += 1) {
      result[i] = this.getLabel(offsets[i]);
    }
    return result;
  }

  abstract size(): number;

  abstract subset(labels: LabelArray): LabelIndexBase;

  abstract isubset(offsets: OffsetArray): LabelIndexBase;

  abstract isubsetMask(mask: Uint8Array | boolean[]): LabelIndexBase;

  abstract withLabel(label: LabelType): LabelIndexBase;

  abstract withLabels(labels: LabelArray): LabelIndexBase;

  abstract dropLabel(label: LabelType): LabelIndexBase;
}

export class IdentityInt32Index extends LabelIndexBase {
  readonly maxOffset: number;

  /*
  identity/noop index, with small assumptions that labels are int32
  */
  constructor(maxOffset: number) {
    super(`IdentityInt32Index_${maxOffset}`);
    this.maxOffset = maxOffset;
  }

  labels(): LabelArray {
    // memoize
    const k = fillRange(new Int32Array(this.maxOffset));
    this.labels = function labels() {
      return k;
    };
    return k;
  }

  getOffset(label: LabelType): number {
    // label to offset
    return Number.isInteger(label) && label >= 0 && label < this.maxOffset
      ? (label as number)
      : -1;
  }

  getLabel(offset: number): number | undefined {
    // offset to label
    return Number.isInteger(offset) && offset >= 0 && offset < this.maxOffset
      ? offset
      : undefined;
  }

  size(): number {
    return this.maxOffset;
  }

  /** @internal */
  __promote(labelArray: LabelArray, allInts: boolean): LabelIndexBase {
    /*
    time/space decision - based on the resulting density
    */
    if (labelArray.length === 0) return new KeyIndex(Array.from(labelArray));
    if (allInts) {
      const [minLabel, maxLabel] = extent(
        labelArray as GenericLabelArray<number> // safe, as allInts is true
      );
      if (minLabel === 0 && maxLabel === labelArray.length - 1)
        return new IdentityInt32Index(labelArray.length);

      const labelSpaceSize = maxLabel - minLabel + 1;
      const density = labelSpaceSize / this.maxOffset;
      /* 0.1 is a magic number which needs testing to optimize */
      if (density < 0.1) {
        return new KeyIndex(Array.from(labelArray));
      }
      return new DenseInt32Index(labelArray as GenericLabelArray<number>, [
        minLabel,
        maxLabel,
      ]);
    }
    return new KeyIndex(Array.from(labelArray));
  }

  subset(labels: LabelArray): LabelIndexBase {
    /* validate subset */
    const { maxOffset } = this;
    for (let i = 0, l = labels.length; i < l; i += 1) {
      const label = labels[i];
      if (!Number.isInteger(label) || label < 0 || label >= maxOffset)
        throw new RangeError(`label: ${label}`);
    }
    return this.__promote(labels, true);
  }

  /* identity index - labels are offsets */
  isubset(offsets: OffsetArray): LabelIndexBase {
    /* validate isubset */
    const { maxOffset } = this;
    for (let i = 0, l = offsets.length; i < l; i += 1) {
      const offset = offsets[i];
      if (!Number.isInteger(offset) || offset < 0 || offset >= maxOffset)
        throw new RangeError(`offset: ${offset}`);
    }
    if (!(offsets instanceof Int32Array)) {
      offsets = new Int32Array(offsets);
    }
    return this.__promote(offsets, true);
  }

  /* identity index - labels are offsets */
  isubsetMask(mask: Uint8Array | boolean[]): LabelIndexBase {
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

  withLabel(label: LabelType): LabelIndexBase {
    if (label === this.maxOffset) {
      return new IdentityInt32Index(label + 1);
    }
    return this.__promote([...this.labels(), label], Number.isInteger(label));
  }

  withLabels(labels: LabelArray): LabelIndexBase {
    return this.__promote(
      [...this.labels(), ...labels],
      labels.every(Number.isInteger)
    );
  }

  dropLabel(label: LabelType): LabelIndexBase {
    if (!Number.isInteger(label) || label < 0 || label > this.maxOffset - 1)
      throw new RangeError("Invalid label.");
    if (label === this.maxOffset - 1) {
      return new IdentityInt32Index(label);
    }
    const labelArray = [...this.labels()];
    labelArray.splice(labelArray.indexOf(label as number), 1);
    return this.__promote(labelArray, true);
  }
}

export class DenseInt32Index extends LabelIndexBase {
  getLabel: (offset: number) => number | undefined;

  getOffset: (label: LabelType) => number;

  index: Int32Array;

  minLabel: number;

  rindex: Int32Array;

  /*
  DenseInt32Index indexes integer labels, and uses Int32Array typed arrays
  for both forward and reverse indexing.   This means that the min/max range
  of the forward index labels must be known a priori (so that the index
  array can be pre-allocated).
  */
  constructor(
    labels: GenericLabelArray<number>,
    labelRange?: [number, number]
  ) {
    super(__getMemoId());
    const int32Labels =
      labels instanceof Int32Array ? labels : new Int32Array(labels);
    if (!labelRange) {
      labelRange = extent(int32Labels);
    }
    const [minLabel, maxLabel] = labelRange;
    const labelSpaceSize = maxLabel - minLabel + 1;
    const index = new Int32Array(labelSpaceSize).fill(-1);
    for (let i = 0, l = labels.length; i < l; i += 1) {
      const label = labels[i];
      index[label - minLabel] = i;
    }

    this.minLabel = minLabel;
    this.rindex = int32Labels;
    this.index = index;

    this.getOffset = function getOffset(label: LabelType) {
      if (!Number.isInteger(label)) return -1;
      const lblIdx: number = <number>label - minLabel;
      if (lblIdx < 0 || lblIdx >= index.length) return -1;
      const offset = index[lblIdx];
      return offset;
    };

    this.getLabel = function getLabel(offset: number) {
      return Number.isInteger(offset) ? labels[offset] : undefined;
    };
  }

  labels(): LabelArray {
    return this.rindex;
  }

  size(): number {
    return this.rindex.length;
  }

  /** @internal */
  __promote(labelArray: LabelArray, allInts: boolean): LabelIndexBase {
    /*
    time/space decision - if we are going to use less than 10% of the
    dense index space, switch to a KeyIndex (which is slower, but uses
    less memory for sparse label spaces).
    */
    if (labelArray.length === 0) return new KeyIndex(Array.from(labelArray));
    if (allInts) {
      if (!(labelArray instanceof Int32Array)) {
        labelArray = new Int32Array(labelArray as number[]);
      }
      const [minLabel, maxLabel] = extent(
        labelArray as GenericLabelArray<number> // safe, as allInts is true
      );
      const labelSpaceSize = maxLabel - minLabel + 1;
      const density = labelSpaceSize / this.rindex.length;
      /* 0.1 is a magic number, that needs testing to optimize */
      if (density < 0.1) {
        return new KeyIndex(Array.from(labelArray));
      }
      return new DenseInt32Index(labelArray as GenericLabelArray<number>, [
        minLabel,
        maxLabel,
      ]);
    }
    return new KeyIndex(Array.from(labelArray));
  }

  subset(labels: LabelArray): LabelIndexBase {
    /* validate subset */
    for (let i = 0, l = labels.length; i < l; i += 1) {
      const label = labels[i]; // if not a number, getOffset will error
      const offset = this.getOffset(label as number);
      if (offset === -1) throw new RangeError(`unknown label: ${label}`);
    }
    return this.__promote(labels as GenericLabelArray<number>, true);
  }

  isubset(offsets: OffsetArray): LabelIndexBase {
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
    return this.__promote(labels, true);
  }

  isubsetMask(mask: Uint8Array | boolean[]): LabelIndexBase {
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
    return this.__promote(labels, true);
  }

  withLabel(label: LabelType): LabelIndexBase {
    return this.__promote([...this.labels(), label], Number.isInteger(label));
  }

  withLabels(labels: LabelArray): LabelIndexBase {
    return this.__promote(
      [...this.labels(), ...labels],
      labels.every(Number.isInteger)
    );
  }

  dropLabel(label: LabelType): LabelIndexBase {
    if (!Number.isInteger(label)) throw new RangeError("Invalid label.");
    const labelArray = [...this.labels()];
    labelArray.splice(labelArray.indexOf(label as number), 1);
    return this.__promote(
      new Int32Array(labelArray as GenericLabelArray<number>),
      true
    );
  }
}

export class KeyIndex extends LabelIndexBase {
  getLabel: (offset: number) => LabelType | undefined;

  getOffset: (label: LabelType) => number | -1;

  index: Map<string | number, number>;

  rindex: (string | number)[];

  /*
  KeyIndex indexes arbitrary JS primitive types, and uses a Map()
  as its core data structure.
  */
  constructor(labels: Array<string | number>) {
    super(__getMemoId());
    const index = new Map<string | number, number>();
    if (labels === undefined) {
      labels = [];
    }
    if (!Array.isArray(labels)) {
      labels = Array.from(labels);
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

    this.getOffset = function getOffset(label: LabelType) {
      const offset = index.get(label);
      if (offset === undefined) return -1;
      return offset;
    };

    this.getLabel = function getLabel(offset: number) {
      return Number.isInteger(offset) ? rindex[offset] : undefined;
    };
  }

  labels(): LabelArray {
    return this.rindex;
  }

  size(): number {
    return this.rindex.length;
  }

  subset(labels: (string | number)[]): LabelIndexBase {
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

  isubset(offsets: OffsetArray): LabelIndexBase {
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

  isubsetMask(mask: Uint8Array | boolean[]): LabelIndexBase {
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

  withLabel(label: LabelType): LabelIndexBase {
    return new KeyIndex([...this.rindex, label]);
  }

  withLabels(labels: LabelArray): LabelIndexBase {
    return new KeyIndex([...this.rindex, ...labels]);
  }

  dropLabel(label: LabelType): LabelIndexBase {
    const idx = this.rindex.indexOf(label);
    const labelArray = [...this.rindex];
    labelArray.splice(idx, 1);
    return new KeyIndex(labelArray);
  }
}

export type LabelIndex = LabelIndexBase;

export function isLabelIndex(i: unknown): i is LabelIndex {
  return (
    i instanceof LabelIndexBase ||
    i instanceof IdentityInt32Index ||
    i instanceof DenseInt32Index ||
    i instanceof KeyIndex
  );
}

/** @internal */
function extent(tarr: GenericLabelArray<number>): [number, number] {
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

/* eslint-enable max-classes-per-file -- enable*/
