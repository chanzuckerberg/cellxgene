// jshint esversion: 6

const BitArray = require("../../src/util/typedCrossfilter/bitArray");
const defaultTestLength = 8;

describe("default select state", () => {
  test("newly created Bitarray should be deselected", () => {
    const ba = new BitArray(defaultTestLength);
    expect(ba).toBeDefined();

    for (let i = 0; i < defaultTestLength; i++) {
      expect(ba.isSelected(i)).toEqual(false);
    }

    const dim = ba.allocDimension();
    expect(dim).toBeDefined();
    for (let i = 0; i < defaultTestLength; i++) {
      expect(ba.isSelected(i)).toEqual(false);
    }

    ba.freeDimension(dim);
    for (let i = 0; i < defaultTestLength; i++) {
      expect(ba.isSelected(i)).toEqual(false);
    }
  });
});

describe("select and deselect", () => {
  test("selectAll and deselectAll", () => {
    const ba = new BitArray(defaultTestLength);
    expect(ba).toBeDefined();
    const dim1 = ba.allocDimension();
    expect(dim1).toBeDefined();
    ba.selectAll(dim1);

    for (let i = 0; i < defaultTestLength; i++) {
      expect(ba.isSelected(i)).toEqual(true);
    }

    const dim2 = ba.allocDimension();
    expect(dim2).toBeDefined();
    for (let i = 0; i < defaultTestLength; i++) {
      expect(ba.isSelected(i)).toEqual(false);
    }

    ba.selectAll(dim2);
    for (let i = 0; i < defaultTestLength; i++) {
      expect(ba.isSelected(i)).toEqual(true);
    }

    ba.deselectAll(dim1);
    for (let i = 0; i < defaultTestLength; i++) {
      expect(ba.isSelected(i)).toEqual(false);
    }

    ba.deselectAll(dim2);
    for (let i = 0; i < defaultTestLength; i++) {
      expect(ba.isSelected(i)).toEqual(false);
    }

    ba.selectAll(dim1);
    ba.selectAll(dim2);
    for (let i = 0; i < defaultTestLength; i++) {
      expect(ba.isSelected(i)).toEqual(true);
    }

    ba.freeDimension(dim1);
    for (let i = 0; i < defaultTestLength; i++) {
      expect(ba.isSelected(i)).toEqual(true);
    }

    ba.freeDimension(dim2);
    for (let i = 0; i < defaultTestLength; i++) {
      expect(ba.isSelected(i)).toEqual(false);
    }
  });

  test("selectOne and deselectOne", () => {
    const ba = new BitArray(defaultTestLength);
    expect(ba).toBeDefined();
    const dim = ba.allocDimension();
    expect(dim).toBeDefined();

    ba.selectOne(dim, 0);
    expect(ba.isSelected(0)).toEqual(true);
    for (let i = 1; i < defaultTestLength; i++) {
      expect(ba.isSelected(i)).toEqual(false);
    }

    ba.deselectOne(dim, 0);
    for (let i = 0; i < defaultTestLength; i++) {
      expect(ba.isSelected(i)).toEqual(false);
    }

    ba.selectOne(dim, 1);
    expect(ba.isSelected(1)).toEqual(true);
    expect(ba.isSelected(0)).toEqual(false);
    for (let i = 2; i < defaultTestLength; i++) {
      expect(ba.isSelected(i)).toEqual(false);
    }

    ba.selectAll(dim);
    console.log(ba);
    console.log(defaultTestLength - 1);
    ba.deselectOne(dim, defaultTestLength - 1);
    console.log(ba);
    expect(ba.isSelected(defaultTestLength - 1)).toEqual(false);
    for (let i = 0; i < defaultTestLength - 1; i++) {
      expect(ba.isSelected(i)).toEqual(true);
    }
  });
});

describe("selectionCount", () => {
  test("simple", () => {
    const ba = new BitArray(defaultTestLength);
    expect(ba).toBeDefined();
    const dim1 = ba.allocDimension();
    expect(dim1).toBeDefined();
    const dim2 = ba.allocDimension();
    expect(dim2).toBeDefined();

    expect(ba.selectionCount).toEqual(0);
    ba.selectAll(dim1);
    expect(ba.selectionCount).toEqual(0);
    ba.selectAll(dim2);
    expect(ba.selectionCount).toEqual(defaultTestLength);

    for (let i = 0; i < defaultTestLength; i++) {
      ba.deselectOne(dim1, i);
      expect(ba.selectionCount).toEqual(defaultTestLength - i - 1);
    }

    ba.freeDimension(dim1);
    ba.freeDimension(dim2);
  });
});
