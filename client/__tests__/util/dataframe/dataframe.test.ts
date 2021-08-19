import * as Dataframe from "../../../src/util/dataframe";

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
      // @ts-expect-error ts-migrate(2345) FIXME: Argument of type 'DenseInt32Index' is not assignab... Remove this comment to see the full error message
      new Dataframe.DenseInt32Index([2, 1, 0]),
      new Dataframe.KeyIndex(["A", "B"])
    );

    expect(df).toBeDefined();
    expect(df.dims).toEqual([3, 2]);

    expect(df.rowIndex).toBeInstanceOf(Dataframe.DenseInt32Index);
    expect(df.colIndex).toBeInstanceOf(Dataframe.KeyIndex);
    expect(df.rowIndex.labels()).toEqual(new Int32Array([2, 1, 0]));
    expect(df.colIndex.labels()).toEqual(["A", "B"]);

    expect(df.at(0, "A")).toEqual(2);
    expect(df.at(2, "B")).toEqual(3);
    expect(df.iat(0, 0)).toEqual(0);
    expect(df.iat(2, 1)).toEqual(5);
  });
});

describe("simple data access", () => {
  const df = new Dataframe.Dataframe(
    [4, 2],
    [
      new Float64Array([0.0, Number.NaN, Number.POSITIVE_INFINITY, 3.14159]),
      ["red", "blue", "green", "nan"],
    ],
    // @ts-expect-error ts-migrate(2345) FIXME: Argument of type 'DenseInt32Index' is not assignab... Remove this comment to see the full error message
    new Dataframe.DenseInt32Index([3, 2, 1, 0]),
    new Dataframe.KeyIndex(["numbers", "colors"])
  );

  test("iat", () => {
    expect(df).toBeDefined();

    // present
    expect(df.iat(0, 0)).toEqual(0.0);
    expect(df.iat(0, 1)).toEqual("red");
    expect(df.iat(1, 0)).toEqual(Number.NaN);
    expect(df.iat(1, 1)).toEqual("blue");
    expect(df.iat(2, 0)).toEqual(Number.POSITIVE_INFINITY);
    expect(df.iat(2, 1)).toEqual("green");
    expect(df.iat(3, 0)).toEqual(3.14159);
    expect(df.iat(3, 1)).toEqual("nan");

    // labels out of range have no defined behavior
  });

  test("at", () => {
    expect(df).toBeDefined();

    // present
    expect(df.at(3, "numbers")).toEqual(0.0);
    expect(df.at(3, "colors")).toEqual("red");
    expect(df.at(2, "numbers")).toEqual(Number.NaN);
    expect(df.at(2, "colors")).toEqual("blue");
    expect(df.at(1, "numbers")).toEqual(Number.POSITIVE_INFINITY);
    expect(df.at(1, "colors")).toEqual("green");
    expect(df.at(0, "numbers")).toEqual(3.14159);
    expect(df.at(0, "colors")).toEqual("nan");

    // labels out of range have no defined behavior
  });

  test("ihas", () => {
    expect(df).toBeDefined();

    // present
    expect(df.ihas(0, 0)).toBeTruthy();
    expect(df.ihas(1, 1)).toBeTruthy();
    expect(df.ihas(3, 1)).toBeTruthy();

    // not present
    expect(df.ihas(-1, -1)).toBeFalsy();
    expect(df.ihas(0, 99)).toBeFalsy();
    expect(df.ihas(99, 0)).toBeFalsy();
    expect(df.ihas(99, 99)).toBeFalsy();
    expect(df.ihas(-1, 0)).toBeFalsy();
    expect(df.ihas(0, -1)).toBeFalsy();
  });

  test("has", () => {
    expect(df).toBeDefined();

    // present
    expect(df.has(3, "numbers")).toBeTruthy();
    expect(df.has(0, "numbers")).toBeTruthy();
    expect(df.has(3, "colors")).toBeTruthy();
    expect(df.has(0, "colors")).toBeTruthy();

    // not present
    expect(df.has(3, "foo")).toBeFalsy();
    expect(df.has(-1, "numbers")).toBeFalsy();
    expect(df.has(-1, -1)).toBeFalsy();
    expect(df.has(null, null)).toBeFalsy();
    expect(df.has(0, "foo")).toBeFalsy();
    expect(df.has(99, "numbers")).toBeFalsy();
    expect(df.has(99, "foo")).toBeFalsy();
  });
});

describe("dataframe subsetting", () => {
  describe("subset", () => {
    const sourceDf = new Dataframe.Dataframe(
      [3, 4],
      [
        new Int32Array([0, 1, 2]),
        ["A", "B", "C"],
        new Float32Array([4.4, 5.5, 6.6]),
        ["red", "green", "blue"],
      ],
      null, // identity index
      // @ts-expect-error ts-migrate(2345) FIXME: Argument of type 'KeyIndex' is not assignable to p... Remove this comment to see the full error message
      new Dataframe.KeyIndex(["int32", "string", "float32", "colors"])
    );

    test("all rows, one column", () => {
      // @ts-expect-error ts-migrate(2345) FIXME: Argument of type 'string[]' is not assignable to p... Remove this comment to see the full error message
      const dfA = sourceDf.subset(null, ["colors"]);
      expect(dfA).toBeDefined();
      expect(dfA.dims).toEqual([3, 1]);
      expect(dfA.iat(0, 0)).toEqual("red");
      expect(dfA.at(2, "colors")).toEqual("blue");
      expect(dfA.col("colors").asArray()).toEqual(["red", "green", "blue"]);
      expect(dfA.icol(0).asArray()).toEqual(["red", "green", "blue"]);
      expect(dfA.col("colors").asArray()).toEqual(
        sourceDf.col("colors").asArray()
      );
      expect(dfA.rowIndex.labels()).toEqual(sourceDf.rowIndex.labels());
      expect(dfA.colIndex.labels()).toEqual(["colors"]);
    });

    test("all rows, two columns", () => {
      // @ts-expect-error ts-migrate(2345) FIXME: Argument of type 'string[]' is not assignable to p... Remove this comment to see the full error message
      const dfB = sourceDf.subset(null, ["float32", "colors"]);
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
      expect(dfB.rowIndex.labels()).toEqual(sourceDf.rowIndex.labels());
      expect(dfB.colIndex.labels()).toEqual(["float32", "colors"]);
    });

    test("one row, all columns", () => {
      const dfC = sourceDf.subset([1], null);
      expect(dfC).toBeDefined();
      expect(dfC.dims).toEqual([1, 4]);
      expect(dfC.iat(0, 0)).toEqual(1);
      expect(dfC.iat(0, 1)).toEqual("B");
      expect(dfC.iat(0, 2)).toBeCloseTo(5.5);
      expect(dfC.iat(0, 3)).toEqual("green");
      expect(dfC.rowIndex.labels()).toEqual(new Int32Array([1]));
      expect(dfC.colIndex.labels()).toEqual(sourceDf.colIndex.labels());
    });

    test("two rows, all columns", () => {
      const dfD = sourceDf.subset([0, 2], null);
      expect(dfD).toBeDefined();
      expect(dfD.dims).toEqual([2, 4]);
      expect(dfD.icol(0).asArray()).toEqual(new Int32Array([0, 2]));
      expect(dfD.icol(1).asArray()).toEqual(["A", "C"]);
      expect(dfD.icol(2).asArray()).toEqual(new Float32Array([4.4, 6.6]));
      expect(dfD.icol(3).asArray()).toEqual(["red", "blue"]);
      expect(dfD.rowIndex.labels()).toEqual(new Int32Array([0, 2]));
      expect(dfD.colIndex.labels()).toEqual(sourceDf.colIndex.labels());

      // reverse the row order
      const dfDr = sourceDf.subset([2, 0], null);
      expect(dfDr.icol(0).asArray()).toEqual(new Int32Array([2, 0]));
      expect(dfDr.icol(1).asArray()).toEqual(["C", "A"]);
      expect(dfDr.icol(2).asArray()).toEqual(new Float32Array([6.6, 4.4]));
      expect(dfDr.icol(3).asArray()).toEqual(["blue", "red"]);
      expect(dfDr.rowIndex.labels()).toEqual(new Int32Array([2, 0]));
      expect(dfDr.colIndex.labels()).toEqual(sourceDf.colIndex.labels());
    });

    test("all rows, all columns", () => {
      const dfE = sourceDf.subset(null, null);
      expect(dfE).toBeDefined();
      expect(dfE.dims).toEqual([3, 4]);
      expect(dfE.icol(0).asArray()).toEqual(sourceDf.icol(0).asArray());
      expect(dfE.icol(1).asArray()).toEqual(sourceDf.icol(1).asArray());
      expect(dfE.icol(2).asArray()).toEqual(sourceDf.icol(2).asArray());
      expect(dfE.icol(3).asArray()).toEqual(sourceDf.icol(3).asArray());
      expect(dfE.rowIndex.labels()).toEqual(sourceDf.rowIndex.labels());
      expect(dfE.colIndex.labels()).toEqual(sourceDf.colIndex.labels());
    });

    test("two rows, two colums", () => {
      // @ts-expect-error ts-migrate(2345) FIXME: Argument of type 'string[]' is not assignable to p... Remove this comment to see the full error message
      const dfF = sourceDf.subset([0, 2], ["int32", "float32"]);
      expect(dfF).toBeDefined();
      expect(dfF.dims).toEqual([2, 2]);
      expect(dfF.icol(0).asArray()).toEqual(new Int32Array([0, 2]));
      expect(dfF.icol(1).asArray()).toEqual(new Float32Array([4.4, 6.6]));
      expect(dfF.rowIndex.labels()).toEqual(new Int32Array([0, 2]));
      expect(dfF.colIndex.labels()).toEqual(["int32", "float32"]);

      // reverse the row and column order
      // @ts-expect-error ts-migrate(2345) FIXME: Argument of type 'string[]' is not assignable to p... Remove this comment to see the full error message
      const dfFr = sourceDf.subset([2, 0], ["float32", "int32"]);
      expect(dfFr).toBeDefined();
      expect(dfFr.dims).toEqual([2, 2]);
      expect(dfFr.icol(0).asArray()).toEqual(new Float32Array([6.6, 4.4]));
      expect(dfFr.icol(1).asArray()).toEqual(new Int32Array([2, 0]));
      expect(dfFr.rowIndex.labels()).toEqual(new Int32Array([2, 0]));
      expect(dfFr.colIndex.labels()).toEqual(["float32", "int32"]);
    });

    test("withRowIndex", () => {
      const df = sourceDf.subset(
        null,
        // @ts-expect-error ts-migrate(2345) FIXME: Argument of type 'string[]' is not assignable to p... Remove this comment to see the full error message
        ["int32", "float32"],
        new Dataframe.DenseInt32Index([3, 2, 1])
      );
      expect(df.colIndex).toBeInstanceOf(Dataframe.KeyIndex);
      expect(df.rowIndex).toBeInstanceOf(Dataframe.DenseInt32Index);
      expect(df.at(3, "int32")).toEqual(df.iat(0, 0));
    });

    test("withRowIndex error checks", () => {
      expect(() =>
        // @ts-expect-error ts-migrate(2345) FIXME: Argument of type 'string[]' is not assignable to p... Remove this comment to see the full error message
        sourceDf.subset(null, ["red"], new Dataframe.IdentityInt32Index(1))
      ).toThrow(RangeError);
      expect(() =>
        // @ts-expect-error ts-migrate(2345) FIXME: Argument of type 'string[]' is not assignable to p... Remove this comment to see the full error message
        sourceDf.subset(null, ["red"], new Dataframe.DenseInt32Index([0, 1]))
      ).toThrow(RangeError);
      expect(() =>
        // @ts-expect-error ts-migrate(2345) FIXME: Argument of type 'string[]' is not assignable to p... Remove this comment to see the full error message
        sourceDf.subset(null, ["red"], new Dataframe.KeyIndex([0, 1, 2, 3]))
      ).toThrow(RangeError);
    });
  });

  test("isubsetMask", () => {
    const sourceDf = new Dataframe.Dataframe(
      [3, 4],
      [
        new Int32Array([0, 1, 2]),
        ["A", "B", "C"],
        new Float32Array([4.4, 5.5, 6.6]),
        ["red", "green", "blue"],
      ],
      // @ts-expect-error ts-migrate(2345) FIXME: Argument of type 'DenseInt32Index' is not assignab... Remove this comment to see the full error message
      new Dataframe.DenseInt32Index([2, 4, 6]),
      new Dataframe.KeyIndex(["int32", "string", "float32", "colors"])
    );

    const dfA = sourceDf.isubsetMask(
      new Uint8Array([0, 1, 1]),
      // @ts-expect-error ts-migrate(2345) FIXME: Argument of type 'Uint8Array' is not assignable to... Remove this comment to see the full error message
      new Uint8Array([1, 0, 0, 1])
    );
    expect(dfA.dims).toEqual([2, 2]);
    expect(dfA.icol(0).asArray()).toEqual(new Int32Array([1, 2]));
    expect(dfA.icol(1).asArray()).toEqual(["green", "blue"]);
    expect(dfA.rowIndex.labels()).toEqual(new Int32Array([4, 6]));
    expect(dfA.colIndex.labels()).toEqual(["int32", "colors"]);
  });

  describe("isubset", () => {
    const sourceDf = new Dataframe.Dataframe(
      [3, 4],
      [
        new Int32Array([0, 1, 2]),
        ["A", "B", "C"],
        new Float32Array([4.4, 5.5, 6.6]),
        ["red", "green", "blue"],
      ],
      null, // identity index
      // @ts-expect-error ts-migrate(2345) FIXME: Argument of type 'KeyIndex' is not assignable to p... Remove this comment to see the full error message
      new Dataframe.KeyIndex(["int32", "string", "float32", "colors"])
    );

    test("one row, all cols", () => {
      const dfA = sourceDf.isubset([1], null);
      expect(dfA.dims).toEqual([1, 4]);
      expect(dfA.icol(0).asArray()).toEqual(new Int32Array([1]));
      expect(dfA.icol(1).asArray()).toEqual(["B"]);
      expect(dfA.icol(2).asArray()).toEqual(new Float32Array([5.5]));
      expect(dfA.icol(3).asArray()).toEqual(["green"]);
    });

    test("all rows, two cols", () => {
      // @ts-expect-error ts-migrate(2345) FIXME: Argument of type 'number[]' is not assignable to p... Remove this comment to see the full error message
      const dfA = sourceDf.isubset(null, [1, 2]);
      expect(dfA.dims).toEqual([3, 2]);
      expect(dfA.icol(0).asArray()).toEqual(["A", "B", "C"]);
      expect(dfA.icol(1).asArray()).toEqual(new Float32Array([4.4, 5.5, 6.6]));
      expect(dfA.col("string")).toBe(dfA.icol(0));
      expect(dfA.col("float32")).toBe(dfA.icol(1));
    });

    test("out of order rows", () => {
      const dfA = sourceDf.isubset([2, 0], null);
      expect(dfA.dims).toEqual([2, 4]);
      expect(dfA.icol(0).asArray()).toEqual(new Int32Array([2, 0]));
      expect(dfA.icol(1).asArray()).toEqual(["C", "A"]);
      expect(dfA.icol(2).asArray()).toEqual(new Float32Array([6.6, 4.4]));
      expect(dfA.icol(3).asArray()).toEqual(["blue", "red"]);
    });
  });
});

describe("dataframe factories", () => {
  test("create", () => {
    const df = Dataframe.Dataframe.create(
      [3, 3],
      [
        new Array(3).fill(0),
        new Int16Array(3).fill(99),
        new Float64Array(3).fill(1.1),
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
      // @ts-expect-error ts-migrate(2345) FIXME: Argument of type 'DenseInt32Index' is not assignab... Remove this comment to see the full error message
      new Dataframe.DenseInt32Index([2, 1, 0]),
      new Dataframe.KeyIndex(["A", "B"])
    );

    const dfB = dfA.clone();
    expect(dfB).not.toBe(dfA);
    expect(dfB.dims).toEqual(dfA.dims);
    expect(dfB).toHaveLength(dfA.length);
    expect(dfB.rowIndex.labels()).toEqual(dfA.rowIndex.labels());
    expect(dfB.colIndex.labels()).toEqual(dfA.colIndex.labels());
    for (let i = 0, l = dfB.dims[1]; i < l; i += 1) {
      expect(dfB.icol(i).asArray()).toEqual(dfA.icol(i).asArray());
    }
  });

  describe("withCol", () => {
    test("KeyIndex", () => {
      const df = new Dataframe.Dataframe(
        [2, 2],
        [
          ["red", "blue"],
          [true, false],
        ],
        null,
        // @ts-expect-error ts-migrate(2345) FIXME: Argument of type 'KeyIndex' is not assignable to p... Remove this comment to see the full error message
        new Dataframe.KeyIndex(["colors", "bools"])
      );
      const dfA = df.withCol("numbers", [1, 0]);

      expect(dfA).toBeDefined();
      expect(dfA.dims).toEqual([2, 3]);
      expect(dfA.icol(0).asArray()).toEqual(["red", "blue"]);
      expect(dfA.icol(1).asArray()).toEqual([true, false]);
      expect(dfA.icol(2).asArray()).toEqual([1, 0]);
      expect(dfA.col("numbers").asArray()).toEqual([1, 0]);
      expect(dfA.colIndex.labels()).toEqual(["colors", "bools", "numbers"]);
      expect(df.colIndex.labels()).toEqual(["colors", "bools"]);
      expect(df.rowIndex.labels()).toEqual(dfA.rowIndex.labels());
    });

    test("DenseInt32Index", () => {
      const df = new Dataframe.Dataframe(
        [2, 2],
        [
          ["red", "blue"],
          [true, false],
        ],
        null,
        // @ts-expect-error ts-migrate(2345) FIXME: Argument of type 'DenseInt32Index' is not assignab... Remove this comment to see the full error message
        new Dataframe.DenseInt32Index([74, 75])
      );
      const dfA = df.withCol(72, [1, 0]);

      expect(dfA).toBeDefined();
      expect(dfA.dims).toEqual([2, 3]);
      expect(dfA.icol(0).asArray()).toEqual(["red", "blue"]);
      expect(dfA.icol(1).asArray()).toEqual([true, false]);
      expect(dfA.icol(2).asArray()).toEqual([1, 0]);
      expect(dfA.col(74).asArray()).toEqual(["red", "blue"]);
      expect(dfA.col(75).asArray()).toEqual([true, false]);
      expect(dfA.col(72).asArray()).toEqual([1, 0]);
      expect(dfA.colIndex.labels()).toEqual(new Int32Array([74, 75, 72]));
      expect(df.colIndex.labels()).toEqual(new Int32Array([74, 75]));
      expect(df.rowIndex.labels()).toEqual(dfA.rowIndex.labels());
    });

    test("DenseInt32Index promote", () => {
      const df = new Dataframe.Dataframe(
        [2, 2],
        [
          ["red", "blue"],
          [true, false],
        ],
        null,
        // @ts-expect-error ts-migrate(2345) FIXME: Argument of type 'DenseInt32Index' is not assignab... Remove this comment to see the full error message
        new Dataframe.DenseInt32Index([74, 75])
      );
      const dfA = df.withCol(999, [1, 0]);

      expect(dfA).toBeDefined();
      expect(dfA.dims).toEqual([2, 3]);
      expect(dfA.icol(0).asArray()).toEqual(["red", "blue"]);
      expect(dfA.icol(1).asArray()).toEqual([true, false]);
      expect(dfA.icol(2).asArray()).toEqual([1, 0]);
      expect(dfA.col(74).asArray()).toEqual(["red", "blue"]);
      expect(dfA.col(75).asArray()).toEqual([true, false]);
      expect(dfA.col(999).asArray()).toEqual([1, 0]);
      expect(dfA.colIndex.labels()).toEqual(new Int32Array([74, 75, 999]));
      expect(df.colIndex.labels()).toEqual(new Int32Array([74, 75]));
      expect(df.rowIndex.labels()).toEqual(dfA.rowIndex.labels());
    });

    test("IdentityInt32Index with last", () => {
      const df = new Dataframe.Dataframe(
        [2, 2],
        [
          ["red", "blue"],
          [true, false],
        ],
        null,
        null
      );
      const dfA = df.withCol(2, [1, 0]);

      expect(dfA).toBeDefined();
      expect(dfA.dims).toEqual([2, 3]);
      expect(dfA.icol(0).asArray()).toEqual(["red", "blue"]);
      expect(dfA.icol(1).asArray()).toEqual([true, false]);
      expect(dfA.icol(2).asArray()).toEqual([1, 0]);
      expect(dfA.col(0).asArray()).toEqual(["red", "blue"]);
      expect(dfA.col(1).asArray()).toEqual([true, false]);
      expect(dfA.col(2).asArray()).toEqual([1, 0]);
      expect(dfA.colIndex.labels()).toEqual(new Int32Array([0, 1, 2]));
      expect(df.colIndex.labels()).toEqual(new Int32Array([0, 1]));
      expect(df.rowIndex.labels()).toEqual(dfA.rowIndex.labels());
    });

    test("IdentityInt32Index promote", () => {
      const df = new Dataframe.Dataframe(
        [2, 2],
        [
          ["red", "blue"],
          [true, false],
        ],
        null,
        null
      );
      const dfA = df.withCol(99, [1, 0]);

      expect(dfA).toBeDefined();
      expect(dfA.dims).toEqual([2, 3]);
      expect(dfA.icol(0).asArray()).toEqual(["red", "blue"]);
      expect(dfA.icol(1).asArray()).toEqual([true, false]);
      expect(dfA.icol(2).asArray()).toEqual([1, 0]);
      expect(dfA.col(0).asArray()).toEqual(["red", "blue"]);
      expect(dfA.col(1).asArray()).toEqual([true, false]);
      expect(dfA.col(99).asArray()).toEqual([1, 0]);
      expect(dfA.colIndex.labels()).toEqual(new Int32Array([0, 1, 99]));
      expect(df.colIndex.labels()).toEqual(new Int32Array([0, 1]));
      expect(df.rowIndex.labels()).toEqual(dfA.rowIndex.labels());
    });

    describe("handle column dimensions correctly", () => {
      /*
      there are two conditions:
      - empty dataframe - will accept an add of any dimensionality
      - non-empty dataframe - added column must match row-count dimension
      */
      test("empty.withCol", () => {
        const edf = Dataframe.Dataframe.empty();
        const df = edf.withCol("foo", [1, 2, 3]);

        expect(edf).toBeDefined();
        expect(df).toBeDefined();
        expect(edf).not.toEqual(df);
        expect(df.dims).toEqual([3, 1]);
        expect(df.icol(0).asArray()).toEqual([1, 2, 3]);
      });

      test("withCol dimension check", () => {
        const dfA = new Dataframe.Dataframe([1, 1], [["a"]]);
        expect(() => {
          dfA.withCol(1, []);
        }).toThrow(RangeError);
      });
    });
  });

  describe("withColsFrom", () => {
    test("error conditions", () => {
      /*
      make sure we catch common errors:
      - duplicate column names
      - dimensionality difference
      */
      const dfA = new Dataframe.Dataframe(
        [2, 3],
        [
          ["red", "blue"],
          [true, false],
          [1, 0],
        ],
        null,
        // @ts-expect-error ts-migrate(2345) FIXME: Argument of type 'KeyIndex' is not assignable to p... Remove this comment to see the full error message
        new Dataframe.KeyIndex(["colors", "bools", "numbers"])
      );

      /* different dimensionality should throw error */
      const dfB = new Dataframe.Dataframe(
        [3, 1],
        [["red", "blue", "green"]],
        null,
        // @ts-expect-error ts-migrate(2345) FIXME: Argument of type 'KeyIndex' is not assignable to p... Remove this comment to see the full error message
        new Dataframe.KeyIndex(["colorsA"])
      );
      // @ts-expect-error ts-migrate(2554) FIXME: Expected 2 arguments, but got 1.
      expect(() => dfA.withColsFrom(dfB)).toThrow(RangeError);

      /* duplicate labels should throw an error */
      // @ts-expect-error ts-migrate(2554) FIXME: Expected 2 arguments, but got 1.
      expect(() => dfA.withColsFrom(dfA)).toThrow(Error);
    });

    test("simple", () => {
      /* simple test that it works as expected in common case */
      const dfEmpty = Dataframe.Dataframe.empty();
      const dfA = new Dataframe.Dataframe(
        [2, 1],
        [["red", "blue"]],
        null,
        // @ts-expect-error ts-migrate(2345) FIXME: Argument of type 'KeyIndex' is not assignable to p... Remove this comment to see the full error message
        new Dataframe.KeyIndex(["colors"])
      );
      const dfB = new Dataframe.Dataframe(
        [2, 1],
        [[true, false]],
        null,
        // @ts-expect-error ts-migrate(2345) FIXME: Argument of type 'KeyIndex' is not assignable to p... Remove this comment to see the full error message
        new Dataframe.KeyIndex(["bools"])
      );

      // @ts-expect-error ts-migrate(2554) FIXME: Expected 2 arguments, but got 1.
      const dfLikeA = dfEmpty.withColsFrom(dfA);
      expect(dfLikeA).toBeDefined();
      expect(dfLikeA.dims).toEqual(dfA.dims);
      expect(dfLikeA.colIndex.labels()).toEqual(dfA.colIndex.labels());
      expect(dfLikeA.rowIndex).toEqual(dfA.rowIndex);
      expect(dfLikeA.rowIndex.labels()).toEqual(dfA.rowIndex.labels());
      expect(dfLikeA.icol(0).asArray()).toEqual(dfA.icol(0).asArray());

      // @ts-expect-error ts-migrate(2554) FIXME: Expected 2 arguments, but got 1.
      const dfAlsoLikeA = dfA.withColsFrom(dfEmpty);
      expect(dfAlsoLikeA).toBeDefined();
      expect(dfAlsoLikeA.dims).toEqual(dfA.dims);
      expect(dfAlsoLikeA.colIndex.labels()).toEqual(dfA.colIndex.labels());
      expect(dfAlsoLikeA.rowIndex).toEqual(dfA.rowIndex);
      expect(dfAlsoLikeA.rowIndex.labels()).toEqual(dfA.rowIndex.labels());
      expect(dfAlsoLikeA.icol(0).asArray()).toEqual(dfA.icol(0).asArray());

      // @ts-expect-error ts-migrate(2554) FIXME: Expected 2 arguments, but got 1.
      const dfC = dfA.withColsFrom(dfB);
      expect(dfC).toBeDefined();
      expect(dfC.dims).toEqual([2, 2]);
      expect(dfC.colIndex.labels()).toEqual(["colors", "bools"]);
      expect(dfC.rowIndex).toEqual(dfA.rowIndex);
      expect(dfC.rowIndex.labels()).toEqual(dfA.rowIndex.labels());
      expect(dfC.col("colors").asArray()).toEqual(["red", "blue"]);
      expect(dfC.col("bools").asArray()).toEqual([true, false]);
    });

    test("column picking", () => {
      const dfEmpty = Dataframe.Dataframe.empty();
      const dfA = new Dataframe.Dataframe(
        [2, 1],
        [["red", "blue"]],
        null,
        // @ts-expect-error ts-migrate(2345) FIXME: Argument of type 'KeyIndex' is not assignable to p... Remove this comment to see the full error message
        new Dataframe.KeyIndex(["colors"])
      );
      const dfB = new Dataframe.Dataframe(
        [2, 3],
        [
          ["red", "blue"],
          [true, false],
          [1, 0],
        ],
        null,
        // @ts-expect-error ts-migrate(2345) FIXME: Argument of type 'KeyIndex' is not assignable to p... Remove this comment to see the full error message
        new Dataframe.KeyIndex(["colors", "bools", "numbers"])
      );

      const dfX = dfEmpty.withColsFrom(dfB, ["colors", "bools"]);
      expect(dfX).toBeDefined();
      expect(dfX.dims).toEqual([2, 2]);
      expect(dfX.colIndex.labels()).toEqual(["colors", "bools"]);
      expect(dfX.rowIndex).toEqual(dfB.rowIndex);
      expect(dfX.icol(0).asArray()).toEqual(dfB.icol(0).asArray());

      const dfY = dfA.withColsFrom(dfB, ["numbers"]);
      expect(dfY).toBeDefined();
      expect(dfY.dims).toEqual([2, 2]);
      expect(dfY.colIndex.labels()).toEqual(["colors", "numbers"]);
      expect(dfY.rowIndex).toEqual(dfA.rowIndex);
      expect(dfY.icol(0).asArray()).toEqual(dfA.icol(0).asArray());

      const dfZ = dfA.withColsFrom(dfEmpty, []);
      expect(dfZ).toBeDefined();
      expect(dfZ.dims).toEqual(dfA.dims);
      expect(dfZ.colIndex.labels()).toEqual(dfA.colIndex.labels());
      expect(dfZ.rowIndex).toEqual(dfA.rowIndex);
      expect(dfZ.icol(0).asArray()).toEqual(dfA.icol(0).asArray());

      expect(() => dfA.withColsFrom(dfB, ["bools", "colors"])).toThrow();
    });

    test("column aliasing", () => {
      const dfA = new Dataframe.Dataframe(
        [2, 1],
        [["red", "blue"]],
        null,
        // @ts-expect-error ts-migrate(2345) FIXME: Argument of type 'KeyIndex' is not assignable to p... Remove this comment to see the full error message
        new Dataframe.KeyIndex(["colors"])
      );
      const dfB = new Dataframe.Dataframe(
        [2, 3],
        [
          ["red", "blue"],
          [true, false],
          [1, 0],
        ],
        null,
        // @ts-expect-error ts-migrate(2345) FIXME: Argument of type 'KeyIndex' is not assignable to p... Remove this comment to see the full error message
        new Dataframe.KeyIndex(["colors", "bools", "numbers"])
      );

      const dfX = dfA.withColsFrom(dfB, { colors: "_colors", bools: "_bools" });
      expect(dfX).toBeDefined();
      expect(dfX.dims).toEqual([2, 3]);
      expect(dfX.colIndex.labels()).toEqual(["colors", "_colors", "_bools"]);
      expect(dfX.rowIndex).toEqual(dfA.rowIndex);
      expect(dfX.icol(0).asArray()).toEqual(dfA.icol(0).asArray());
      expect(dfX.col("_colors").asArray()).toBe(dfB.col("colors").asArray());
    });
  });

  describe("dropCol", () => {
    test("KeyIndex", () => {
      const df = new Dataframe.Dataframe(
        [2, 3],
        [
          ["red", "blue"],
          [true, false],
          [1, 0],
        ],
        null,
        // @ts-expect-error ts-migrate(2345) FIXME: Argument of type 'KeyIndex' is not assignable to p... Remove this comment to see the full error message
        new Dataframe.KeyIndex(["colors", "bools", "numbers"])
      );
      const dfA = df.dropCol("colors");

      expect(dfA).toBeDefined();
      expect(dfA.dims).toEqual([2, 2]);
      expect(dfA.icol(0).asArray()).toEqual([true, false]);
      expect(dfA.icol(1).asArray()).toEqual([1, 0]);
      expect(dfA.col("numbers").asArray()).toEqual([1, 0]);
      expect(dfA.colIndex.labels()).toEqual(["bools", "numbers"]);
      expect(df.colIndex.labels()).toEqual(["colors", "bools", "numbers"]);
      expect(df.rowIndex.labels()).toEqual(dfA.rowIndex.labels());
    });

    test("IdentityInt32Index drop first", () => {
      const df = new Dataframe.Dataframe(
        [2, 3],
        [
          ["red", "blue"],
          [true, false],
          [1, 0],
        ],
        null,
        null
      );
      const dfA = df.dropCol(0);

      expect(dfA).toBeDefined();
      expect(dfA.dims).toEqual([2, 2]);
      expect(dfA.icol(0).asArray()).toEqual([true, false]);
      expect(dfA.icol(1).asArray()).toEqual([1, 0]);
      expect(df.col(1).asArray()).toEqual(dfA.col(1).asArray());
      expect(df.col(2).asArray()).toEqual(dfA.col(2).asArray());
      expect(dfA.colIndex.labels()).toEqual(new Int32Array([1, 2]));
      expect(df.colIndex.labels()).toEqual(new Int32Array([0, 1, 2]));
      expect(df.rowIndex.labels()).toEqual(dfA.rowIndex.labels());
    });

    test("IdentityInt32Index drop last", () => {
      const df = new Dataframe.Dataframe(
        [2, 3],
        [
          ["red", "blue"],
          [true, false],
          [1, 0],
        ],
        null,
        null
      );
      const dfA = df.dropCol(2);

      expect(dfA).toBeDefined();
      expect(dfA.dims).toEqual([2, 2]);
      expect(dfA.icol(0).asArray()).toEqual(["red", "blue"]);
      expect(dfA.icol(1).asArray()).toEqual([true, false]);
      expect(df.col(0).asArray()).toEqual(dfA.col(0).asArray());
      expect(df.col(1).asArray()).toEqual(dfA.col(1).asArray());
      expect(dfA.colIndex.labels()).toEqual(new Int32Array([0, 1]));
      expect(df.colIndex.labels()).toEqual(new Int32Array([0, 1, 2]));
      expect(df.rowIndex.labels()).toEqual(dfA.rowIndex.labels());
    });

    test("DenseInt32Index", () => {
      const df = new Dataframe.Dataframe(
        [2, 3],
        [
          ["red", "blue"],
          [true, false],
          [1, 0],
        ],
        null,
        // @ts-expect-error ts-migrate(2345) FIXME: Argument of type 'DenseInt32Index' is not assignab... Remove this comment to see the full error message
        new Dataframe.DenseInt32Index([102, 101, 100])
      );
      const dfA = df.dropCol(101);

      expect(dfA).toBeDefined();
      expect(dfA.dims).toEqual([2, 2]);
      expect(dfA.icol(0).asArray()).toEqual(["red", "blue"]);
      expect(dfA.icol(1).asArray()).toEqual([1, 0]);
      expect(dfA.col(100).asArray()).toEqual([1, 0]);
      expect(dfA.col(102).asArray()).toEqual(["red", "blue"]);
      expect(dfA.colIndex.labels()).toEqual(new Int32Array([102, 100]));
      expect(df.colIndex.labels()).toEqual(new Int32Array([102, 101, 100]));
      expect(df.rowIndex.labels()).toEqual(dfA.rowIndex.labels());
    });
  });

  describe("mapColumns", () => {
    test("identity", () => {
      const dfA = Dataframe.Dataframe.create(
        [3, 3],
        [
          new Array(3).fill(0),
          new Int16Array(3).fill(99),
          new Float64Array(3).fill(1.1),
        ]
      );
      // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
      const dfB = dfA.mapColumns((col: any, idx: any) => {
        expect(dfA.icol(idx).asArray()).toBe(col);
        return col;
      });
      expect(dfA).not.toBe(dfB);
      expect(dfA.dims).toEqual(dfB.dims);
      for (let c = 0; c < dfA.dims[1]; c += 1) {
        expect(dfA.icol(c).asArray()).toBe(dfB.icol(c).asArray());
      }
    });

    test("transform", () => {
      const dfA = Dataframe.Dataframe.create(
        [3, 3],
        [new Array(3).fill(0), new Array(3).fill(0), new Array(3).fill(0)]
      );
      const dfB = dfA.mapColumns(() => new Array(3).fill(1));
      expect(dfA).not.toBe(dfB);
      expect(dfB.iat(0, 0)).toEqual(1);
      expect(dfB.iat(0, 1)).toEqual(1);
      expect(dfB.iat(0, 2)).toEqual(1);
    });

    test("columns", () => {
      const df = Dataframe.Dataframe.create(
        [3, 3],
        [new Array(3).fill(0), new Array(3).fill(0), new Array(3).fill(0)]
      );

      expect(df).toBeDefined();
      expect(df.columns()).toHaveLength(3);
      expect(df.columns()[0]).toEqual(df.icol(0));
      expect(df.columns()[2]).toEqual(df.icol(2));
    });

    test("renameCol", () => {
      const dfA = new Dataframe.Dataframe(
        [2, 2],
        [
          [true, false],
          [1, 0],
        ],
        null,
        // @ts-expect-error ts-migrate(2345) FIXME: Argument of type 'KeyIndex' is not assignable to p... Remove this comment to see the full error message
        new Dataframe.KeyIndex(["A", "B"])
      );
      const dfB = dfA.renameCol("B", "C");
      expect(dfA.colIndex.labels()).toEqual(["A", "B"]);
      expect(dfB.colIndex.labels()).toEqual(["A", "C"]);
      expect(dfA.dims).toMatchObject(dfB.dims);
      expect(dfA.columns()).toMatchObject(dfB.columns());
    });
  });
});

describe("dataframe col", () => {
  // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
  let df: any = null;
  beforeEach(() => {
    df = new Dataframe.Dataframe(
      [2, 2],
      [
        [true, false],
        [1, 0],
      ],
      null,
      // @ts-expect-error ts-migrate(2345) FIXME: Argument of type 'KeyIndex' is not assignable to p... Remove this comment to see the full error message
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

describe("label indexing", () => {
  describe("isLabelIndex", () => {
    expect(
      Dataframe.isLabelIndex(new Dataframe.IdentityInt32Index(4))
    ).toBeTruthy();
    expect(
      Dataframe.isLabelIndex(new Dataframe.DenseInt32Index([2, 4, 99]))
    ).toBeTruthy();
    expect(
      Dataframe.isLabelIndex(new Dataframe.KeyIndex(["a", 4, "toasty"]))
    ).toBeTruthy();
    expect(Dataframe.isLabelIndex(false)).toBeFalsy();
    expect(Dataframe.isLabelIndex(undefined)).toBeFalsy();
    expect(Dataframe.isLabelIndex(null)).toBeFalsy();
    expect(Dataframe.isLabelIndex(true)).toBeFalsy();
    expect(Dataframe.isLabelIndex([])).toBeFalsy();
    expect(Dataframe.isLabelIndex({})).toBeFalsy();
    expect(Dataframe.isLabelIndex(Dataframe.IdentityInt32Index)).toBeFalsy();
  });

  describe("IdentityInt32Index", () => {
    const idx = new Dataframe.IdentityInt32Index(12); // [0, 12)

    test("create", () => {
      expect(Dataframe.isLabelIndex(idx)).toBeTruthy();
    });

    test("labels", () => {
      expect(idx.labels()).toEqual(
        new Int32Array([0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11])
      );
      expect(idx.getLabel(1)).toEqual(1);
      expect(idx.getLabels([1, 3])).toEqual([1, 3]);
      expect(idx.size()).toEqual(12);
    });

    test("offsets", () => {
      expect(idx.getOffset(1)).toEqual(1);
      expect(idx.getOffsets([1, 3])).toEqual([1, 3]);
    });

    test("subset", () => {
      expect(idx.subset([2]).labels()).toEqual([2]);
      expect(idx.subset([2, 3, 4]).labels()).toEqual(new Int32Array([2, 3, 4]));
      expect(idx.subset([0, 1, 2, 3]).labels()).toEqual(
        new Int32Array([0, 1, 2, 3])
      );
      expect(idx.subset([0, 1, 2, 3, 4])).toBeInstanceOf(
        Dataframe.IdentityInt32Index
      );
      expect(idx.subset([2, 1, 0])).toBeInstanceOf(
        Dataframe.IdentityInt32Index
      );
      expect(idx.subset([1, 2, 3, 4])).toBeInstanceOf(
        Dataframe.DenseInt32Index
      );
      expect(idx.subset([0, 1, 3, 4])).toBeInstanceOf(
        Dataframe.DenseInt32Index
      );
      expect(idx.subset([0, 1, 2, 3, 10])).toBeInstanceOf(
        Dataframe.DenseInt32Index
      );
      expect(idx.subset([4, 3, 2, 1])).toBeInstanceOf(
        Dataframe.DenseInt32Index
      );
      expect(idx.subset([4])).toBeInstanceOf(Dataframe.KeyIndex);
    });

    test("isubset", () => {
      expect(idx.isubset([2]).labels()).toEqual([2]);
      expect(idx.isubset([2, 3, 4]).labels()).toEqual(
        new Int32Array([2, 3, 4])
      );
      expect(idx.isubset([0, 1, 2, 3]).labels()).toEqual(
        new Int32Array([0, 1, 2, 3])
      );
      expect(() => idx.isubset([-1001])).toThrow(RangeError);
      expect(() => idx.isubset([1001])).toThrow(RangeError);
    });

    test("isubsetMask", () => {
      expect(
        idx
          .isubsetMask([
            true,
            true,
            true,
            true,
            true,
            true,
            true,
            true,
            true,
            true,
            true,
            true,
          ])
          .labels()
      ).toEqual(new Int32Array([0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11]));
      expect(
        idx
          .isubsetMask([
            false,
            false,
            false,
            false,
            false,
            false,
            false,
            false,
            false,
            false,
            false,
            false,
          ])
          .labels()
      ).toEqual(new Int32Array([]));
      expect(
        idx
          .isubsetMask([
            false,
            true,
            true,
            true,
            true,
            true,
            true,
            true,
            true,
            true,
            true,
            true,
          ])
          .labels()
      ).toEqual(new Int32Array([1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11]));
      expect(
        idx
          .isubsetMask([
            false,
            true,
            true,
            true,
            true,
            true,
            true,
            true,
            true,
            true,
            false,
            true,
          ])
          .labels()
      ).toEqual(new Int32Array([1, 2, 3, 4, 5, 6, 7, 8, 9, 11]));
      expect(
        idx
          .isubsetMask([
            false,
            true,
            true,
            false,
            true,
            true,
            true,
            true,
            true,
            true,
            true,
            false,
          ])
          .labels()
      ).toEqual(new Int32Array([1, 2, 4, 5, 6, 7, 8, 9, 10]));
      expect(() => idx.isubsetMask([])).toThrow(RangeError);
    });

    test("withLabel", () => {
      expect(idx.withLabel(99).labels()).toEqual(
        new Int32Array([0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 99])
      );
      expect(idx.withLabel(12).labels()).toEqual(
        new Int32Array([0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12])
      );
      expect(idx.withLabels([12, 13]).labels()).toEqual(
        new Int32Array([0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13])
      );
    });

    test("dropLabel", () => {
      expect(idx.dropLabel(0).labels()).toEqual(
        new Int32Array([1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11])
      );
      expect(idx.dropLabel(11).labels()).toEqual(
        new Int32Array([0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10])
      );
      expect(idx.dropLabel(5).labels()).toEqual(
        new Int32Array([0, 1, 2, 3, 4, 6, 7, 8, 9, 10, 11])
      );
    });
  });

  describe("DenseInt32Index", () => {
    const idx = new Dataframe.DenseInt32Index([99, 1002, 48, 0, 22]);

    test("create", () => {
      expect(Dataframe.isLabelIndex(idx)).toBeTruthy();
    });

    test("labels", () => {
      expect(idx.labels()).toEqual(new Int32Array([99, 1002, 48, 0, 22]));
      expect(idx.size()).toEqual(5);
      expect(idx.getLabel(0)).toEqual(99);
      expect(idx.getLabels(new Int32Array([2, 4]))).toEqual(
        new Int32Array([48, 22])
      );
      expect(idx.getLabels([2, 4])).toEqual([48, 22]);
    });

    test("offsets", () => {
      expect(idx.getOffset(1002)).toEqual(1);
      expect(idx.getOffset(0)).toEqual(3);
      expect(idx.getOffsets([0, 48])).toEqual([3, 2]);
    });

    test("subset", () => {
      expect(idx.subset([1002, 0, 99]).labels()).toEqual(
        new Int32Array([1002, 0, 99])
      );
      expect(idx.getOffsets(idx.subset([1002, 0, 99]).labels())).toEqual(
        new Int32Array([1, 3, 0])
      );
      expect(() => idx.subset([-1])).toThrow(RangeError);
    });

    test("isubset", () => {
      expect(idx.isubset([4, 1, 2]).labels()).toEqual(
        new Int32Array([22, 1002, 48])
      );
      expect(() => idx.isubset([-1001])).toThrow(RangeError);
      expect(() => idx.isubset([1001])).toThrow(RangeError);
    });

    test("isubsetMask", () => {
      expect(idx.isubsetMask([true, true, true, true, true]).labels()).toEqual(
        new Int32Array([99, 1002, 48, 0, 22])
      );
      expect(
        idx.isubsetMask([false, false, false, false, false]).labels()
      ).toEqual(new Int32Array([]));
      expect(idx.isubsetMask([true, true, false, true, true]).labels()).toEqual(
        new Int32Array([99, 1002, 0, 22])
      );
      expect(
        idx.isubsetMask([false, true, true, true, false]).labels()
      ).toEqual(new Int32Array([1002, 48, 0]));
      expect(() => idx.isubsetMask([])).toThrow(RangeError);
    });

    test("withLabel", () => {
      expect(idx.withLabel(88).labels()).toEqual(
        new Int32Array([99, 1002, 48, 0, 22, 88])
      );
      expect(idx.withLabel(88).getOffset(88)).toEqual(5);
      expect(idx.withLabels([88, 99]).labels()).toEqual(
        new Int32Array([99, 1002, 48, 0, 22, 88, 99])
      );
    });
    test("dropLabel", () => {
      expect(idx.dropLabel(48).labels()).toEqual(
        new Int32Array([99, 1002, 0, 22])
      );
    });
  });

  describe("KeyIndex", () => {
    const idx = new Dataframe.KeyIndex(["red", "green", "blue"]);

    test("create", () => {
      expect(Dataframe.isLabelIndex(idx)).toBeTruthy();
      expect(() => new Dataframe.KeyIndex(["dup", "dup"])).toThrow(Error);
      // @ts-expect-error ts-migrate(2554) FIXME: Expected 1 arguments, but got 0.
      expect(new Dataframe.KeyIndex().size()).toEqual(0);
    });

    test("labels", () => {
      expect(idx.labels()).toEqual(["red", "green", "blue"]);
      expect(idx.size()).toEqual(3);
      expect(idx.getLabel(1)).toEqual("green");
      expect(idx.getLabels([2, 0])).toEqual(["blue", "red"]);
    });

    test("offsets", () => {
      expect(idx.getOffset("blue")).toEqual(2);
    });

    test("subset", () => {
      expect(idx.subset(["green"]).labels()).toEqual(["green"]);
      expect(idx.subset(["green", "red"]).labels()).toEqual(["green", "red"]);
    });

    test("isubset", () => {
      expect(idx.isubset([2, 1, 0]).labels()).toEqual(["blue", "green", "red"]);
      expect(() => idx.isubset([-1001])).toThrow(RangeError);
      expect(() => idx.isubset([1001])).toThrow(RangeError);
    });

    test("isubsetMask", () => {
      expect(idx.isubsetMask([true, true, true]).labels()).toEqual([
        "red",
        "green",
        "blue",
      ]);
      expect(idx.isubsetMask([false, false, false]).labels()).toEqual([]);
      expect(idx.isubsetMask([true, false, true]).labels()).toEqual([
        "red",
        "blue",
      ]);
      expect(() => idx.isubsetMask([])).toThrow(RangeError);
    });

    test("withLabel", () => {
      expect(idx.withLabel("yo").labels()).toEqual([
        "red",
        "green",
        "blue",
        "yo",
      ]);
      expect(idx.withLabel("yo").getOffset("yo")).toEqual(3);
      expect(idx.withLabels(["hey", "there"]).labels()).toEqual([
        "red",
        "green",
        "blue",
        "hey",
        "there",
      ]);
    });

    test("dropLabel", () => {
      expect(idx.dropLabel("blue").labels()).toEqual(["red", "green"]);
    });
  });
});

describe("corner cases", () => {
  /* error/corner cases */

  test("identity integer index rejects non-integer labels", () => {
    const idx = new Dataframe.IdentityInt32Index(10);
    expect(idx.getOffset(0)).toBe(0);
    expect(idx.getOffset(9)).toBe(9);
    expect(idx.getOffset(10)).toBeUndefined();
    expect(idx.getOffset(-1)).toBeUndefined();
    expect(idx.getOffset("sort")).toBeUndefined();
    expect(idx.getOffset("length")).toBeUndefined();
    expect(idx.getOffset(true)).toBeUndefined();
    expect(idx.getOffset(0.001)).toBeUndefined();
    expect(idx.getOffset({})).toBeUndefined();
    expect(idx.getOffset([])).toBeUndefined();
    expect(idx.getOffset(new Float32Array())).toBeUndefined();
    expect(idx.getOffset("__proto__")).toBeUndefined();

    expect(idx.getLabel(0)).toBe(0);
    expect(idx.getLabel(9)).toBe(9);
    expect(idx.getLabel(10)).toBeUndefined();
    expect(idx.getLabel(-1)).toBeUndefined();
    expect(idx.getLabel("sort")).toBeUndefined();
    expect(idx.getLabel("length")).toBeUndefined();
    expect(idx.getLabel(true)).toBeUndefined();
    expect(idx.getLabel(0.001)).toBeUndefined();
    expect(idx.getLabel({})).toBeUndefined();
    expect(idx.getLabel([])).toBeUndefined();
    expect(idx.getLabel(new Float32Array())).toBeUndefined();
    expect(idx.getLabel("__proto__")).toBeUndefined();
  });

  test("dense integer index rejects non-integer labels", () => {
    const idx = new Dataframe.DenseInt32Index([-10, 0, 3, 9, 10]);
    expect(idx.getOffset(0)).toBe(1);
    expect(idx.getOffset(9)).toBe(3);
    expect(idx.getOffset(1)).toBeUndefined();
    expect(idx.getOffset(11)).toBeUndefined();
    expect(idx.getOffset(-1)).toBeUndefined();
    expect(idx.getOffset("sort")).toBeUndefined();
    expect(idx.getOffset("length")).toBeUndefined();
    expect(idx.getOffset(true)).toBeUndefined();
    expect(idx.getOffset(0.001)).toBeUndefined();
    expect(idx.getOffset({})).toBeUndefined();
    expect(idx.getOffset([])).toBeUndefined();
    expect(idx.getOffset(new Float32Array())).toBeUndefined();
    expect(idx.getOffset("__proto__")).toBeUndefined();

    expect(idx.getLabel(0)).toBe(-10);
    expect(idx.getLabel(4)).toBe(10);
    expect(idx.getLabel(10)).toBeUndefined();
    expect(idx.getLabel(-1)).toBeUndefined();
    expect(idx.getLabel("sort")).toBeUndefined();
    expect(idx.getLabel("length")).toBeUndefined();
    expect(idx.getLabel(true)).toBeUndefined();
    expect(idx.getLabel(0.001)).toBeUndefined();
    expect(idx.getLabel({})).toBeUndefined();
    expect(idx.getLabel([])).toBeUndefined();
    expect(idx.getLabel(new Float32Array())).toBeUndefined();
    expect(idx.getLabel("__proto__")).toBeUndefined();
  });

  test("Empty dataframe rejects bogus labels", () => {
    const df = Dataframe.Dataframe.empty();

    expect(df.hasCol("sort")).toBeFalsy();
    expect(df.hasCol(0)).toBeFalsy();
    expect(df.hasCol(true)).toBeFalsy();
    expect(df.hasCol(false)).toBeFalsy();
    expect(df.hasCol([])).toBeFalsy();
    expect(df.hasCol({})).toBeFalsy();
    expect(df.hasCol(null)).toBeFalsy();
    expect(df.hasCol(undefined)).toBeFalsy();

    expect(df.col("sort")).toBeUndefined();
    expect(df.col(0)).toBeUndefined();
    expect(df.col(true)).toBeUndefined();
    expect(df.col(false)).toBeUndefined();
    expect(df.col([])).toBeUndefined();
    expect(df.col({})).toBeUndefined();
    expect(df.col(null)).toBeUndefined();
    expect(df.col(undefined)).toBeUndefined();

    expect(df.icol("sort")).toBeUndefined();
    expect(df.icol(0)).toBeUndefined();
    expect(df.icol(true)).toBeUndefined();
    expect(df.icol(false)).toBeUndefined();
    expect(df.icol([])).toBeUndefined();
    expect(df.icol({})).toBeUndefined();
    expect(df.icol(null)).toBeUndefined();
    expect(df.icol(undefined)).toBeUndefined();

    expect(df.ihas("sort", "length")).toBeFalsy();
    expect(df.ihas("0", "0")).toBeFalsy();
    expect(df.ihas("", "")).toBeFalsy();
    expect(df.ihas(null, null)).toBeFalsy();
    expect(df.ihas(undefined, undefined)).toBeFalsy();
    expect(df.ihas(true, true)).toBeFalsy();
    expect(df.ihas([], [])).toBeFalsy();
    expect(df.ihas({}, {})).toBeFalsy();
  });

  test("Dataframe rejects bogus labels", () => {
    const df = new Dataframe.Dataframe(
      [2, 2],
      [
        [true, false],
        [1, 0],
      ],
      null,
      // @ts-expect-error ts-migrate(2345) FIXME: Argument of type 'KeyIndex' is not assignable to p... Remove this comment to see the full error message
      new Dataframe.KeyIndex(["A", "B"])
    );

    expect(df.hasCol("sort")).toBeFalsy();
    expect(df.hasCol("__proto__")).toBeFalsy();
    expect(df.hasCol(0)).toBeFalsy();
    expect(df.hasCol(true)).toBeFalsy();
    expect(df.hasCol(false)).toBeFalsy();
    expect(df.hasCol([])).toBeFalsy();
    expect(df.hasCol({})).toBeFalsy();
    expect(df.hasCol(null)).toBeFalsy();
    expect(df.hasCol(undefined)).toBeFalsy();

    expect(df.col("sort")).toBeUndefined();
    expect(df.col("__proto__")).toBeUndefined();
    expect(df.col(0)).toBeUndefined();
    expect(df.col(true)).toBeUndefined();
    expect(df.col(false)).toBeUndefined();
    expect(df.col([])).toBeUndefined();
    expect(df.col({})).toBeUndefined();
    expect(df.col(null)).toBeUndefined();
    expect(df.col(undefined)).toBeUndefined();

    expect(df.icol("sort")).toBeUndefined();
    expect(df.icol("__proto__")).toBeUndefined();
    expect(df.icol(-1)).toBeUndefined();
    expect(df.icol(true)).toBeUndefined();
    expect(df.icol(false)).toBeUndefined();
    expect(df.icol([])).toBeUndefined();
    expect(df.icol({})).toBeUndefined();
    expect(df.icol(null)).toBeUndefined();
    expect(df.icol(undefined)).toBeUndefined();

    expect(df.ihas("sort", "length")).toBeFalsy();
    expect(df.ihas("__proto__", "__proto__")).toBeFalsy();

    expect(df.ihas(-1, 0)).toBeFalsy();
    expect(df.ihas("0", 0)).toBeFalsy();
    expect(df.ihas("", 0)).toBeFalsy();
    expect(df.ihas(null, 0)).toBeFalsy();
    expect(df.ihas(undefined, 0)).toBeFalsy();
    expect(df.ihas([], 0)).toBeFalsy();
    expect(df.ihas({}, 0)).toBeFalsy();

    expect(df.ihas(0, -1)).toBeFalsy();
    expect(df.ihas(0, "0")).toBeFalsy();
    expect(df.ihas(0, "")).toBeFalsy();
    expect(df.ihas(0, null)).toBeFalsy();
    expect(df.ihas(0, undefined)).toBeFalsy();
    expect(df.ihas(0, [])).toBeFalsy();
    expect(df.ihas(0, {})).toBeFalsy();

    expect(df.has("sort", "length")).toBeFalsy();
    expect(df.has("length", "sort")).toBeFalsy();
    expect(df.has("__proto__", "__proto__")).toBeFalsy();

    expect(df.has(-1, "A")).toBeFalsy();
    expect(df.has("0", "A")).toBeFalsy();
    expect(df.has("", "A")).toBeFalsy();
    expect(df.has(null, "A")).toBeFalsy();
    expect(df.has(undefined, "A")).toBeFalsy();
    expect(df.has([], "A")).toBeFalsy();
    expect(df.has({}, "A")).toBeFalsy();

    expect(df.has(0, -1)).toBeFalsy();
    expect(df.has(0, "0")).toBeFalsy();
    expect(df.has(0, "")).toBeFalsy();
    expect(df.has(0, null)).toBeFalsy();
    expect(df.has(0, undefined)).toBeFalsy();
    expect(df.has(0, [])).toBeFalsy();
    expect(df.has(0, {})).toBeFalsy();
  });
});
