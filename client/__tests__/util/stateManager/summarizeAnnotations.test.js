import summarizeAnnotations from "../../../src/util/stateManager/summarizeAnnotations";

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
    const summary = summarizeAnnotations(schema, [], []);
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
              max: Number.NEGATIVE_INFINITY,
              min: Number.POSITIVE_INFINITY
            }
          },
          nameInt32: {
            categorical: false,
            range: {
              max: Number.NEGATIVE_INFINITY,
              min: Number.POSITIVE_INFINITY
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
    const obsAnnotations = [
      {
        __index__: 0,
        name: "n1",
        nameString: "hi",
        nameBoolean: true,
        nameFloat32: 39.3,
        nameInt32: 99,
        nameCategorical: 1
      }
    ];
    const varAnnotations = [];

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
            range: { min: 39.3, max: 39.3 }
          },
          nameInt32: {
            categorical: false,
            range: { min: 99, max: 99 }
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
    const obsAnnotations = [
      {
        __index__: 0,
        name: "n0",
        nameString: "hi",
        nameBoolean: false,
        nameFloat32: 39.3,
        nameInt32: 99,
        nameCategorical: 1
      },
      {
        __index__: 1,
        name: "n1",
        nameString: "hi",
        nameBoolean: true,
        nameFloat32: 39.3,
        nameInt32: 99,
        nameCategorical: false
      },
      {
        __index__: 2,
        name: "n2",
        nameString: "bye",
        nameBoolean: true,
        nameFloat32: 0,
        nameInt32: 99,
        nameCategorical: "0"
      }
    ];
    const varAnnotations = [];

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
            range: { min: 0, max: 39.3 }
          },
          nameInt32: {
            categorical: false,
            range: { min: 99, max: 99 }
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
