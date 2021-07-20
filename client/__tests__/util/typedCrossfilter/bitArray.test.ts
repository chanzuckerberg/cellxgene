import BitArray from "../../../src/util/typedCrossfilter/bitArray";

const defaultTestLength = 8;

describe("default select state", () => {
  test("newly created Bitarray should be deselected", () => {
    const ba = new BitArray(defaultTestLength);
    expect(ba).toBeDefined();

    for (let i = 0; i < defaultTestLength; i += 1) {
      expect(ba.isSelected(i)).toEqual(false);
    }

    const dim = ba.allocDimension();
    expect(dim).toBeDefined();
    for (let i = 0; i < defaultTestLength; i += 1) {
      expect(ba.isSelected(i)).toEqual(false);
    }

    ba.freeDimension(dim);
    for (let i = 0; i < defaultTestLength; i += 1) {
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

    for (let i = 0; i < defaultTestLength; i += 1) {
      expect(ba.isSelected(i)).toEqual(true);
    }

    const dim2 = ba.allocDimension();
    expect(dim2).toBeDefined();
    for (let i = 0; i < defaultTestLength; i += 1) {
      expect(ba.isSelected(i)).toEqual(false);
    }

    ba.selectAll(dim2);
    for (let i = 0; i < defaultTestLength; i += 1) {
      expect(ba.isSelected(i)).toEqual(true);
    }

    ba.deselectAll(dim1);
    for (let i = 0; i < defaultTestLength; i += 1) {
      expect(ba.isSelected(i)).toEqual(false);
    }

    ba.deselectAll(dim2);
    for (let i = 0; i < defaultTestLength; i += 1) {
      expect(ba.isSelected(i)).toEqual(false);
    }

    ba.selectAll(dim1);
    ba.selectAll(dim2);
    for (let i = 0; i < defaultTestLength; i += 1) {
      expect(ba.isSelected(i)).toEqual(true);
    }

    ba.freeDimension(dim1);
    for (let i = 0; i < defaultTestLength; i += 1) {
      expect(ba.isSelected(i)).toEqual(true);
    }

    ba.freeDimension(dim2);
    for (let i = 0; i < defaultTestLength; i += 1) {
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
    for (let i = 1; i < defaultTestLength; i += 1) {
      expect(ba.isSelected(i)).toEqual(false);
    }

    ba.deselectOne(dim, 0);
    for (let i = 0; i < defaultTestLength; i += 1) {
      expect(ba.isSelected(i)).toEqual(false);
    }

    ba.selectOne(dim, 1);
    expect(ba.isSelected(1)).toEqual(true);
    expect(ba.isSelected(0)).toEqual(false);
    for (let i = 2; i < defaultTestLength; i += 1) {
      expect(ba.isSelected(i)).toEqual(false);
    }

    ba.selectAll(dim);
    ba.deselectOne(dim, defaultTestLength - 1);
    expect(ba.isSelected(defaultTestLength - 1)).toEqual(false);
    for (let i = 0; i < defaultTestLength - 1; i += 1) {
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

    expect(ba.selectionCount()).toEqual(0);
    ba.selectAll(dim1);
    expect(ba.selectionCount()).toEqual(0);
    ba.selectAll(dim2);
    expect(ba.selectionCount()).toEqual(defaultTestLength);

    for (let i = 0; i < defaultTestLength; i += 1) {
      ba.deselectOne(dim1, i);
      expect(ba.selectionCount()).toEqual(defaultTestLength - i - 1);
      expect(ba.selectionCount()).toEqual(ba.countAllOnes());
    }

    ba.freeDimension(dim1);
    ba.freeDimension(dim2);
  });
});

describe("fillBySelection", () => {
  test("sets values correctly", () => {
    const ba = new BitArray(defaultTestLength);
    expect(ba).toBeDefined();
    const dim = ba.allocDimension();
    expect(dim).toBeDefined();

    const arr = new Int32Array(defaultTestLength);
    arr.fill(0);
    const truth = new Int32Array(defaultTestLength);
    truth.fill(0);

    // initial state should be deselected
    ba.fillBySelection(arr, 1, 0);
    expect(arr).toEqual(expect.not.arrayContaining([1]));

    // selectAll
    ba.selectAll(dim);
    ba.fillBySelection(arr, 1, 0);
    expect(arr).toEqual(expect.not.arrayContaining([0]));

    // deselectOne
    ba.deselectOne(dim, 3);
    ba.fillBySelection(arr, 1, 0);
    truth.fill(1);
    truth[3] = 0;
    expect(arr).toEqual(truth);

    // deselectAll
    ba.deselectAll(dim);
    ba.fillBySelection(arr, 1, 0);
    truth.fill(0);
    expect(arr).toEqual(truth);

    // selectOne
    ba.selectOne(dim, 5);
    ba.fillBySelection(arr, 6, 1);
    truth.fill(1);
    truth[5] = 6;
    expect(arr).toEqual(truth);

    // should be deselected after dimension disposal
    ba.freeDimension(dim);
    ba.fillBySelection(arr, 3, 9);
    truth.fill(9);
    expect(arr).toEqual(truth);
  });
});

describe("wide bitarray", () => {
  test.each([9, 30, 31, 32, 33, 54, 63, 64, 65, 127, 128, 129])(
    "more than %d dimensions",
    (ndim) => {
      /* ensure we move across the uint boundary correctly */

      const ba = new BitArray(defaultTestLength);
      expect(ba).toBeDefined();

      for (let i = 0; i < ndim; i += 1) {
        expect(ba.allocDimension()).toEqual(i);
      }
      ba.freeDimension(0);
      expect(ba.allocDimension()).toEqual(0);

      ba.freeDimension(ndim - 1);
      expect(ba.allocDimension()).toEqual(ndim - 1);
    }
  );
});
