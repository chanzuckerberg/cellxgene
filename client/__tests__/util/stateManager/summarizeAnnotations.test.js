import summarizeAnnotations from "../../../src/util/stateManager/summarizeAnnotations";
import * as Dataframe from "../../../src/util/dataframe";

function float32Conversion(f) {
  return new Float32Array([39.3])[0];
}

describe("summarizeAnnotations", () => {
  const schema = {
    annotations: {
      obs: [
        { name: "name", type: "string" },
        { name: "nameString", type: "string" },
        { name: "nameBoolean", type: "boolean" },
        { name: "nameFloat32", type: "float32" },
        { name: "nameInt32", type: "int32" },
        {
          name: "nameCategorical",
          type: "categorical",
          categories: [true, false, 1, 0, 0.00001, 4383.4833, "test", "", "0"]
        }
      ],
      var: [{ name: "name", type: "string" }]
    }
  };

  test("empty test", () => {
    const df = Dataframe.Dataframe.empty();
    const summary = summarizeAnnotations(schema, df, df.clone());
    expect(summary).toEqual(
      expect.objectContaining({
        obs: {
          nameString: {
            categorical: true,
            categories: [],
            categoryCounts: new Map(),
            numCategories: 0
          },
          nameBoolean: {
            categorical: true,
            categories: [],
            categoryCounts: new Map(),
            numCategories: 0
          },
          nameFloat32: {
            categorical: false,
            range: {
              max: undefined,
              min: undefined,
              nan: 0,
              ninf: 0,
              pinf: 0
            }
          },
          nameInt32: {
            categorical: false,
            range: {
              max: undefined,
              min: undefined,
              nan: 0,
              ninf: 0,
              pinf: 0
            }
          },
          nameCategorical: {
            categorical: true,
            categories: [],
            categoryCounts: new Map(),
            numCategories: 0
          }
        },
        var: {}
      })
    );
  });

  test("simple test", () => {
    const obsAnnotations = new Dataframe.Dataframe(
      [1, 6],
      [
        ["n1"],
        ["hi"],
        [true],
        new Float32Array([39.3]),
        new Int32Array([99]),
        [1]
      ],
      null,
      new Dataframe.KeyIndex([
        "name",
        "nameString",
        "nameBoolean",
        "nameFloat32",
        "nameInt32",
        "nameCategorical"
      ])
    );
    const varAnnotations = Dataframe.Dataframe.empty();

    const summary = summarizeAnnotations(
      schema,
      obsAnnotations,
      varAnnotations
    );

    expect(summary).toEqual(
      expect.objectContaining({
        obs: {
          nameString: {
            categorical: true,
            categories: ["hi"],
            categoryCounts: new Map([["hi", 1]]),
            numCategories: 1
          },
          nameBoolean: {
            categorical: true,
            categories: [true],
            categoryCounts: new Map([[true, 1]]),
            numCategories: 1
          },
          nameFloat32: {
            categorical: false,
            range: {
              min: float32Conversion(39.3),
              max: float32Conversion(39.3),
              nan: 0,
              ninf: 0,
              pinf: 0
            }
          },
          nameInt32: {
            categorical: false,
            range: { min: 99, max: 99, nan: 0, ninf: 0, pinf: 0 }
          },
          nameCategorical: {
            categorical: true,
            categories: [1],
            categoryCounts: new Map([[1, 1]]),
            numCategories: 1
          }
        },
        var: {}
      })
    );
  });

  test("multi test", () => {
    const obsAnnotations = new Dataframe.Dataframe(
      [3, 6],
      [
        ["n0", "n1", "n2"],
        ["hi", "hi", "bye"],
        [false, true, true],
        new Float32Array([39.3, 39.3, 0]),
        new Int32Array([99, 99, 99]),
        [1, false, "0"]
      ],
      null,
      new Dataframe.KeyIndex([
        "name",
        "nameString",
        "nameBoolean",
        "nameFloat32",
        "nameInt32",
        "nameCategorical"
      ])
    );
    const varAnnotations = Dataframe.Dataframe.empty();

    const summary = summarizeAnnotations(
      schema,
      obsAnnotations,
      varAnnotations
    );

    expect(summary).toMatchObject(
      expect.objectContaining({
        obs: {
          nameString: {
            categorical: true,
            categories: expect.arrayContaining(["hi", "bye"]),
            categoryCounts: new Map([["hi", 2], ["bye", 1]]),
            numCategories: 2
          },
          nameBoolean: {
            categorical: true,
            categories: expect.arrayContaining([true, false]),
            categoryCounts: new Map([[true, 2], [false, 1]]),
            numCategories: 2
          },
          nameFloat32: {
            categorical: false,
            range: {
              min: 0,
              max: float32Conversion(39.3),
              nan: 0,
              ninf: 0,
              pinf: 0
            }
          },
          nameInt32: {
            categorical: false,
            range: { min: 99, max: 99, nan: 0, ninf: 0, pinf: 0 }
          },
          nameCategorical: {
            categorical: true,
            categories: expect.arrayContaining([1, false, "0"]),
            categoryCounts: new Map([[1, 1], [false, 1], ["0", 1]]),
            numCategories: 3
          }
        },
        var: {}
      })
    );
  });

  test("non-finite numbers", () => {
    const obsAnnotations = new Dataframe.Dataframe(
      [4, 6],
      [
        ["n0", "n1", "n2", "n2"],
        ["hi", "hi", "bye", "bye"],
        [false, true, true, true],
        new Float32Array([
          39.3,
          Number.NEGATIVE_INFINITY,
          Number.NaN,
          Number.POSITIVE_INFINITY
        ]),
        new Int32Array([99, 99, 99, 99]),
        [1, false, "0", "0"]
      ],
      null,
      new Dataframe.KeyIndex([
        "name",
        "nameString",
        "nameBoolean",
        "nameFloat32",
        "nameInt32",
        "nameCategorical"
      ])
    );
    const varAnnotations = Dataframe.Dataframe.empty();

    const summary = summarizeAnnotations(
      schema,
      obsAnnotations,
      varAnnotations
    );

    expect(summary).toMatchObject(
      expect.objectContaining({
        obs: {
          nameString: {
            categorical: true,
            categories: expect.arrayContaining(["hi", "bye"]),
            categoryCounts: new Map([["hi", 2], ["bye", 1]]),
            numCategories: 2
          },
          nameBoolean: {
            categorical: true,
            categories: expect.arrayContaining([true, false]),
            categoryCounts: new Map([[true, 2], [false, 1]]),
            numCategories: 2
          },
          nameFloat32: {
            categorical: false,
            range: {
              min: float32Conversion(39.3),
              max: float32Conversion(39.3),
              nan: 1,
              ninf: 1,
              pinf: 1
            }
          },
          nameInt32: {
            categorical: false,
            range: { min: 99, max: 99, nan: 0, ninf: 0, pinf: 0 }
          },
          nameCategorical: {
            categorical: true,
            categories: expect.arrayContaining([1, false, "0"]),
            categoryCounts: new Map([[1, 1], [false, 1], ["0", 1]]),
            numCategories: 3
          }
        },
        var: {}
      })
    );
  });
});
