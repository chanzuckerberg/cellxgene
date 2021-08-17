import * as Dataframe from "../../../src/util/dataframe";

// eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
function float32Conversion(f: any) {
  return new Float32Array([f])[0];
}

describe("Dataframe column summary", () => {
  test("empty column test", () => {
    const df = Dataframe.Dataframe.create([0, 1], [[]]);
    const summary = df.icol(0).summarizeCategorical();
    expect(summary).toEqual(
      expect.objectContaining({
        categorical: true,
        categories: [],
        categoryCounts: new Map(),
        numCategories: 0,
      })
    );
  });

  test("simple test", () => {
    const df = new Dataframe.Dataframe(
      [1, 6],
      [
        ["n1"],
        ["hi"],
        [true],
        new Float32Array([39.3]),
        new Int32Array([99]),
        [1],
      ],
      null,
      new Dataframe.KeyIndex([
        "name",
        "nameString",
        "nameBoolean",
        "nameFloat32",
        "nameInt32",
        "nameCategorical",
      ])
    );

    expect(df.icol(0).summarizeCategorical()).toEqual(
      expect.objectContaining({
        categorical: true,
        categories: ["n1"],
        categoryCounts: new Map([["n1", 1]]),
        numCategories: 1,
      })
    );
    expect(df.icol(1).summarizeCategorical()).toEqual(
      expect.objectContaining({
        categorical: true,
        categories: ["hi"],
        categoryCounts: new Map([["hi", 1]]),
        numCategories: 1,
      })
    );
    expect(df.icol(2).summarizeCategorical()).toEqual(
      expect.objectContaining({
        categorical: true,
        categories: [true],
        categoryCounts: new Map([[true, 1]]),
        numCategories: 1,
      })
    );
    expect(df.icol(3).summarizeContinuous()).toEqual(
      expect.objectContaining({
        categorical: false,
        min: float32Conversion(39.3),
        max: float32Conversion(39.3),
        nan: 0,
        ninf: 0,
        pinf: 0,
      })
    );
    expect(df.icol(4).summarizeContinuous()).toEqual(
      expect.objectContaining({
        categorical: false,
        min: 99,
        max: 99,
        nan: 0,
        ninf: 0,
        pinf: 0,
      })
    );
    expect(df.icol(5).summarizeCategorical()).toEqual(
      expect.objectContaining({
        categorical: true,
        categories: [1],
        categoryCounts: new Map([[1, 1]]),
        numCategories: 1,
      })
    );
  });

  test("multi test", () => {
    const df = new Dataframe.Dataframe(
      [3, 6],
      [
        ["n0", "n1", "n2"],
        ["hi", "hi", "bye"],
        [false, true, true],
        new Float32Array([39.3, 39.3, 0]),
        new Int32Array([99, 99, 99]),
        [1, false, "0"],
      ],
      null,
      new Dataframe.KeyIndex([
        "name",
        "nameString",
        "nameBoolean",
        "nameFloat32",
        "nameInt32",
        "nameCategorical",
      ])
    );

    expect(df.icol(0).summarizeCategorical()).toEqual(
      expect.objectContaining({
        categorical: true,
        categories: expect.arrayContaining(["n0", "n1", "n2"]),
        categoryCounts: new Map([
          ["n0", 1],
          ["n1", 1],
          ["n2", 1],
        ]),
        numCategories: 3,
      })
    );
    expect(df.icol(1).summarizeCategorical()).toEqual(
      expect.objectContaining({
        categorical: true,
        categories: expect.arrayContaining(["hi", "bye"]),
        categoryCounts: new Map([
          ["hi", 2],
          ["bye", 1],
        ]),
        numCategories: 2,
      })
    );
    expect(df.icol(2).summarizeCategorical()).toEqual(
      expect.objectContaining({
        categorical: true,
        categories: expect.arrayContaining([true, false]),
        categoryCounts: new Map([
          [true, 2],
          [false, 1],
        ]),
        numCategories: 2,
      })
    );
    expect(df.icol(3).summarizeContinuous()).toEqual(
      expect.objectContaining({
        categorical: false,
        min: 0,
        max: float32Conversion(39.3),
        nan: 0,
        ninf: 0,
        pinf: 0,
      })
    );
    expect(df.icol(4).summarizeContinuous()).toEqual(
      expect.objectContaining({
        categorical: false,
        min: 99,
        max: 99,
        nan: 0,
        ninf: 0,
        pinf: 0,
      })
    );
    expect(df.icol(5).summarizeCategorical()).toEqual(
      expect.objectContaining({
        categorical: true,
        categories: expect.arrayContaining([1, false, "0"]),
        // @ts-expect-error ts-migrate(2769) FIXME: No overload matches this call.
        categoryCounts: new Map([
          [1, 1],
          [false, 1],
          ["0", 1],
        ]),
        numCategories: 3,
      })
    );
  });

  test("non-finite numbers", () => {
    const df = new Dataframe.Dataframe(
      [4, 6],
      [
        ["n0", "n1", "n2", "n2"],
        ["hi", "hi", "bye", "bye"],
        [false, true, true, true],
        new Float32Array([
          39.3,
          Number.NEGATIVE_INFINITY,
          Number.NaN,
          Number.POSITIVE_INFINITY,
        ]),
        new Int32Array([99, 99, 99, 99]),
        [1, false, "0", "0"],
      ],
      null,
      new Dataframe.KeyIndex([
        "name",
        "nameString",
        "nameBoolean",
        "nameFloat32",
        "nameInt32",
        "nameCategorical",
      ])
    );

    expect(df.icol(0).summarizeCategorical()).toEqual(
      expect.objectContaining({
        categorical: true,
        categories: expect.arrayContaining(["n0", "n1", "n2"]),
        categoryCounts: new Map([
          ["n0", 1],
          ["n1", 1],
          ["n2", 2],
        ]),
        numCategories: 3,
      })
    );
    expect(df.icol(1).summarizeCategorical()).toEqual(
      expect.objectContaining({
        categorical: true,
        categories: expect.arrayContaining(["hi", "bye"]),
        categoryCounts: new Map([
          ["hi", 2],
          ["bye", 1],
        ]),
        numCategories: 2,
      })
    );
    expect(df.icol(2).summarizeCategorical()).toEqual(
      expect.objectContaining({
        categorical: true,
        categories: expect.arrayContaining([true, false]),
        categoryCounts: new Map([
          [true, 2],
          [false, 1],
        ]),
        numCategories: 2,
      })
    );
    expect(df.icol(3).summarizeContinuous()).toEqual(
      expect.objectContaining({
        categorical: false,
        min: float32Conversion(39.3),
        max: float32Conversion(39.3),
        nan: 1,
        ninf: 1,
        pinf: 1,
      })
    );
    expect(df.icol(4).summarizeContinuous()).toEqual(
      expect.objectContaining({
        categorical: false,
        min: 99,
        max: 99,
        nan: 0,
        ninf: 0,
        pinf: 0,
      })
    );
    expect(df.icol(5).summarizeCategorical()).toEqual(
      expect.objectContaining({
        categorical: true,
        categories: expect.arrayContaining([1, false, "0"]),
        // @ts-expect-error ts-migrate(2769) FIXME: No overload matches this call.
        categoryCounts: new Map([
          [1, 1],
          [false, 1],
          ["0", 1],
        ]),
        numCategories: 3,
      })
    );
  });
});
