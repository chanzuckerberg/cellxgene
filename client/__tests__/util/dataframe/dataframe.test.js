import * as Dataframe from "../../../src/util/dataframe";

/*
Known gaps in these tests:
- row indexing, especially a la 'world', are not sufficiently tested
*/

describe("dataframe constructor", () => {
  test("empty dataframe", () => {
    const df = new Dataframe.Dataframe([0, 0], []);
    expect(df).toBeDefined();
    expect(df.dims).toEqual([0, 0]);
    expect(df).toHaveLength(0);
    expect(df.icol(0)).not.toBeDefined();
  });

  test("create with default indices", () => {
    const df = new Dataframe.Dataframe(
      [3, 2],
      [new Int32Array(3).fill(0), new Int32Array(3).fill(1)]
    );

    expect(df).toBeDefined();
    expect(df.dims).toEqual([3, 2]);
    expect(df.rowIndex).toBeInstanceOf(Dataframe.IdentityInt32Index);
    expect(df.colIndex).toBeInstanceOf(Dataframe.IdentityInt32Index);
    expect(df.at(0, 0)).toEqual(0);
    expect(df.at(2, 1)).toEqual(1);
    expect(df.iat(0, 0)).toEqual(0);
    expect(df.iat(2, 1)).toEqual(1);
  });

  test("create with labelled indices", () => {
    const df = new Dataframe.Dataframe(
      [3, 2],
      [new Int32Array([0, 1, 2]), new Int32Array([3, 4, 5])],
      new Dataframe.DenseInt32Index([2, 1, 0]),
      new Dataframe.KeyIndex(["A", "B"])
    );

    expect(df).toBeDefined();
    expect(df.dims).toEqual([3, 2]);

    expect(df.rowIndex).toBeInstanceOf(Dataframe.DenseInt32Index);
    expect(df.colIndex).toBeInstanceOf(Dataframe.KeyIndex);
    expect(df.rowIndex.keys()).toEqual(new Int32Array([2, 1, 0]));
    expect(df.colIndex.keys()).toEqual(["A", "B"]);

    expect(df.at(0, "A")).toEqual(2);
    expect(df.at(2, "B")).toEqual(3);
    expect(df.iat(0, 0)).toEqual(0);
    expect(df.iat(2, 1)).toEqual(5);
  });
});

describe("dataframe subsetting", () => {
  test("cutByList", () => {
    const sourceDf = new Dataframe.Dataframe(
      [3, 4],
      [
        new Int32Array([0, 1, 2]),
        ["A", "B", "C"],
        new Float32Array([4.4, 5.5, 6.6]),
        ["red", "green", "blue"]
      ],
      null,
      new Dataframe.KeyIndex(["int32", "string", "float32", "colors"])
    );

    // all rows, one column
    const dfA = sourceDf.cutByList(null, ["colors"]);
    expect(dfA).toBeDefined();
    expect(dfA.dims).toEqual([3, 1]);
    expect(dfA.iat(0, 0)).toEqual("red");
    expect(dfA.at(2, "colors")).toEqual("blue");
    expect(dfA.col("colors").asArray()).toEqual(["red", "green", "blue"]);
    expect(dfA.icol(0).asArray()).toEqual(["red", "green", "blue"]);
    expect(dfA.col("colors").asArray()).toEqual(
      sourceDf.col("colors").asArray()
    );
    expect(dfA.rowIndex.keys()).toEqual(sourceDf.rowIndex.keys());
    expect(dfA.colIndex.keys()).toEqual(["colors"]);

    // all rows, two columns
    const dfB = sourceDf.cutByList(null, ["colors", "float32"]);
    expect(dfB).toBeDefined();
    expect(dfB.dims).toEqual([3, 2]);
    expect(dfB.iat(0, 0)).toBeCloseTo(4.4);
    expect(dfB.iat(0, 1)).toEqual("red");
    expect(dfB.at(2, "colors")).toEqual("blue");
    expect(dfB.at(2, "float32")).toBeCloseTo(6.6);
    expect(dfB.col("colors").asArray()).toEqual(["red", "green", "blue"]);
    expect(dfB.col("float32").asArray()).toEqual(
      new Float32Array([4.4, 5.5, 6.6])
    );
    expect(dfB.icol(0).asArray()).toEqual(dfB.col("float32").asArray());
    expect(dfB.icol(1).asArray()).toEqual(dfB.col("colors").asArray());
    expect(dfB.col("colors").asArray()).toEqual(
      sourceDf.col("colors").asArray()
    );
    expect(dfB.col("float32").asArray()).toEqual(
      sourceDf.col("float32").asArray()
    );
    expect(dfB.rowIndex.keys()).toEqual(sourceDf.rowIndex.keys());
    expect(dfB.colIndex.keys()).toEqual(["float32", "colors"]);

    // one row, all columns
    const dfC = sourceDf.cutByList([1], null);
    expect(dfC).toBeDefined();
    expect(dfC.dims).toEqual([1, 4]);
    expect(dfC.iat(0, 0)).toEqual(1);
    expect(dfC.iat(0, 1)).toEqual("B");
    expect(dfC.iat(0, 2)).toBeCloseTo(5.5);
    expect(dfC.iat(0, 3)).toEqual("green");
    expect(dfC.rowIndex.keys()).toEqual(new Int32Array([1]));
    expect(dfC.colIndex.keys()).toEqual(sourceDf.colIndex.keys());

    // two rows, all columns
    const dfD = sourceDf.cutByList([0, 2], null);
    expect(dfD).toBeDefined();
    expect(dfD.dims).toEqual([2, 4]);
    expect(dfD.icol(0).asArray()).toEqual(new Int32Array([0, 2]));
    expect(dfD.icol(1).asArray()).toEqual(["A", "C"]);
    expect(dfD.icol(2).asArray()).toEqual(new Float32Array([4.4, 6.6]));
    expect(dfD.icol(3).asArray()).toEqual(["red", "blue"]);
    expect(dfD.rowIndex.keys()).toEqual(new Int32Array([0, 2]));
    expect(dfD.colIndex.keys()).toEqual(sourceDf.colIndex.keys());

    // all rows, all columns
    const dfE = sourceDf.cutByList(null, null);
    expect(dfE).toBeDefined();
    expect(dfE.dims).toEqual([3, 4]);
    expect(dfE.icol(0).asArray()).toEqual(sourceDf.icol(0).asArray());
    expect(dfE.icol(1).asArray()).toEqual(sourceDf.icol(1).asArray());
    expect(dfE.icol(2).asArray()).toEqual(sourceDf.icol(2).asArray());
    expect(dfE.icol(3).asArray()).toEqual(sourceDf.icol(3).asArray());
    expect(dfE.rowIndex.keys()).toEqual(sourceDf.rowIndex.keys());
    expect(dfE.colIndex.keys()).toEqual(sourceDf.colIndex.keys());

    // two rows, two colums
    const dfF = sourceDf.cutByList([0, 2], ["int32", "float32"]);
    expect(dfF).toBeDefined();
    expect(dfF.dims).toEqual([2, 2]);
    expect(dfF.icol(0).asArray()).toEqual(new Int32Array([0, 2]));
    expect(dfF.icol(1).asArray()).toEqual(new Float32Array([4.4, 6.6]));
    expect(dfF.rowIndex.keys()).toEqual(new Int32Array([0, 2]));
    expect(dfF.colIndex.keys()).toEqual(["int32", "float32"]);
  });

  test("icutByMask", () => {
    const sourceDf = new Dataframe.Dataframe(
      [3, 4],
      [
        new Int32Array([0, 1, 2]),
        ["A", "B", "C"],
        new Float32Array([4.4, 5.5, 6.6]),
        ["red", "green", "blue"]
      ],
      new Dataframe.DenseInt32Index([2, 4, 6]),
      new Dataframe.KeyIndex(["int32", "string", "float32", "colors"])
    );

    const dfA = sourceDf.icutByMask(
      new Uint8Array([0, 1, 1]),
      new Uint8Array([1, 0, 0, 1])
    );
    expect(dfA.dims).toEqual([2, 2]);
    expect(dfA.icol(0).asArray()).toEqual(new Int32Array([1, 2]));
    expect(dfA.icol(1).asArray()).toEqual(["green", "blue"]);
    expect(dfA.rowIndex.keys()).toEqual(new Int32Array([4, 6]));
    expect(dfA.colIndex.keys()).toEqual(["int32", "colors"]);
  });
});

describe("dataframe factories", () => {
  test("create", () => {
    const df = Dataframe.Dataframe.create(
      [3, 3],
      [
        new Array(3).fill(0),
        new Int16Array(3).fill(99),
        new Float64Array(3).fill(1.1)
      ]
    );

    expect(df).toBeDefined();
    expect(df.dims).toEqual([3, 3]);
    expect(df).toHaveLength(3);
    expect(df.iat(0, 0)).toEqual(0);
    expect(df.iat(1, 1)).toEqual(99);
    expect(df.iat(2, 2)).toBeCloseTo(1.1);
    expect(df.iat(0, 0)).toEqual(df.at(0, 0));
    expect(df.iat(1, 1)).toEqual(df.at(1, 1));
    expect(df.iat(2, 2)).toEqual(df.at(2, 2));
  });

  test("clone", () => {
    const dfA = new Dataframe.Dataframe(
      [3, 2],
      [new Int32Array([0, 1, 2]), new Int32Array([3, 4, 5])],
      new Dataframe.DenseInt32Index([2, 1, 0]),
      new Dataframe.KeyIndex(["A", "B"])
    );

    const dfB = dfA.clone();
    expect(dfB).not.toBe(dfA);
    expect(dfB.dims).toEqual(dfA.dims);
    expect(dfB).toHaveLength(dfA.length);
    expect(dfB.rowIndex.keys()).toEqual(dfA.rowIndex.keys());
    expect(dfB.colIndex.keys()).toEqual(dfA.colIndex.keys());
    for (let i = 0, l = dfB.dims[1]; i < l; i += 1) {
      expect(dfB.icol(i).asArray()).toEqual(dfA.icol(i).asArray());
    }
  });
});

describe("dataframe col", () => {
  let df = null;
  beforeEach(() => {
    df = new Dataframe.Dataframe(
      [2, 2],
      [[true, false], [1, 0]],
      null,
      new Dataframe.KeyIndex(["A", "B"])
    );
  });

  test("col", () => {
    expect(df).toBeDefined();
    expect(df.col("A")).toBe(df.icol(0));
    expect(df.col("B")).toBe(df.icol(1));
    expect(df.col("undefined")).toBeUndefined();
    expect(df.icol("undefined")).toBeUndefined();

    const colA = df.col("A");
    expect(colA).toBeInstanceOf(Function);
    expect(colA.asArray).toBeInstanceOf(Function);
    expect(colA.has).toBeInstanceOf(Function);
    expect(colA.ihas).toBeInstanceOf(Function);
    expect(colA.indexOf).toBeInstanceOf(Function);
    expect(colA.iget).toBeInstanceOf(Function);
  });

  test("col.asArray", () => {
    expect(df).toBeDefined();
    expect(df.col("A").asArray()).toEqual([true, false]);
    expect(df.icol(0).asArray()).toEqual([true, false]);
    expect(df.col("B").asArray()).toEqual([1, 0]);
    expect(df.icol(1).asArray()).toEqual([1, 0]);
  });

  test("col.has", () => {
    expect(df).toBeDefined();
    expect(df.col("A").has(-1)).toBe(false);
    expect(df.col("A").has(0)).toBe(true);
    expect(df.col("A").has(1)).toBe(true);
    expect(df.col("A").has(2)).toBe(false);
    expect(df.col("B").has(-1)).toBe(false);
    expect(df.col("B").has(0)).toBe(true);
    expect(df.col("B").has(1)).toBe(true);
    expect(df.col("B").has(2)).toBe(false);
  });

  test("col.ihas", () => {
    expect(df).toBeDefined();
    expect(df.col("A").ihas(-1)).toBe(false);
    expect(df.col("A").ihas(0)).toBe(true);
    expect(df.col("A").ihas(1)).toBe(true);
    expect(df.col("A").ihas(2)).toBe(false);
    expect(df.col("B").ihas(-1)).toBe(false);
    expect(df.col("B").ihas(0)).toBe(true);
    expect(df.col("B").ihas(1)).toBe(true);
    expect(df.col("B").ihas(2)).toBe(false);
  });

  test("col.iget", () => {
    expect(df).toBeDefined();
    expect(df.col("A").iget(0)).toEqual(df.iat(0, 0));
    expect(df.col("B").iget(1)).toEqual(df.iat(1, 1));
  });

  test("col.indexOf", () => {
    expect(df).toBeDefined();
    expect(df.col("A").indexOf(true)).toEqual(0);
    expect(df.col("A").indexOf(false)).toEqual(1);
    expect(df.col("A").indexOf(99)).toBeUndefined();
    expect(df.col("A").indexOf(undefined)).toBeUndefined();
    expect(df.col("A").indexOf(1)).toBeUndefined();

    expect(df.col("B").indexOf(1)).toEqual(0);
    expect(df.col("B").indexOf(0)).toEqual(1);
    expect(df.col("B").indexOf(99)).toBeUndefined();
    expect(df.col("B").indexOf(undefined)).toBeUndefined();
    expect(df.col("B").indexOf(true)).toBeUndefined();
  });
});
