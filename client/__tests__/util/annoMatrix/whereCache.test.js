import sha1 from "sha1";
import {
  _whereCacheGet,
  _whereCacheCreate,
  _whereCacheMerge,
} from "../../../src/annoMatrix/whereCache";

const schema = {};

describe("whereCache", () => {
  test("whereCacheGet - where query, missing cache values", () => {
    expect(
      _whereCacheGet({}, schema, "X", {
        where: {
          field: "var",
          column: "foo",
          value: "bar",
        },
      })
    ).toEqual([undefined]);
    expect(
      _whereCacheGet({}, schema, "X", {
        summarize: {
          field: "var",
          column: "foo",
          values: ["bar"],
        },
      })
    ).toEqual([undefined]);
    expect(
      _whereCacheGet({ where: { X: {} } }, schema, "X", {
        where: {
          field: "var",
          column: "foo",
          value: "bar",
        },
      })
    ).toEqual([undefined]);
    expect(
      _whereCacheGet({ where: { X: { var: new Map() } } }, schema, "X", {
        where: {
          field: "var",
          column: "foo",
          value: "bar",
        },
      })
    ).toEqual([undefined]);
    expect(
      _whereCacheGet(
        { where: { X: { var: new Map([["foo", new Map()]]) } } },
        schema,
        "X",
        {
          where: {
            field: "var",
            column: "foo",
            value: "bar",
          },
        }
      )
    ).toEqual([undefined]);
  });

  test("whereCacheGet - summarize query, missing cache values", () => {
    expect(
      _whereCacheGet({}, schema, "X", {
        summarize: {
          method: "mean",
          field: "var",
          column: "foo",
          values: ["bar"],
        },
      })
    ).toEqual([undefined]);
    expect(
      _whereCacheGet(
        { summarize: { X: { mean: { var: new Map() } } } },
        schema,
        "X",
        {
          summarize: {
            method: "mean",
            field: "var",
            column: "foo",
            values: ["bar"],
          },
        }
      )
    ).toEqual([undefined]);
  });

  test("whereCacheGet - varied lookups", () => {
    const whereCache = {
      where: {
        X: {
          var: new Map([
            [
              "foo",
              new Map([
                ["bar", [0]],
                ["baz", [1, 2]],
              ]),
            ],
          ]),
        },
      },
      summarize: {
        X: {
          mean: {
            var: new Map([
              [
                "foo",
                new Map([
                  [sha1("bar"), [0]],
                  [sha1("baz"), [1, 2]],
                ]),
              ],
            ]),
          },
        },
      },
    };

    expect(
      _whereCacheGet(whereCache, schema, "X", {
        where: {
          field: "var",
          column: "foo",
          value: "bar",
        },
      })
    ).toEqual([0]);
    expect(
      _whereCacheGet(whereCache, schema, "X", {
        summarize: {
          method: "mean",
          field: "var",
          column: "foo",
          values: ["bar"],
        },
      })
    ).toEqual([0]);
    expect(
      _whereCacheGet(whereCache, schema, "X", {
        where: {
          field: "var",
          column: "foo",
          value: "baz",
        },
      })
    ).toEqual([1, 2]);
    expect(
      _whereCacheGet(whereCache, schema, "X", {
        summarize: {
          method: "mean",
          field: "var",
          column: "foo",
          values: ["baz"],
        },
      })
    ).toEqual([1, 2]);
    expect(_whereCacheGet(whereCache, schema, "Y", {})).toEqual([undefined]);
    expect(
      _whereCacheGet(whereCache, schema, "X", {
        where: {
          field: "whoknows",
          column: "whatever",
          value: "snork",
        },
      })
    ).toEqual([undefined]);
    expect(
      _whereCacheGet(whereCache, schema, "X", {
        where: {
          field: "var",
          column: "whatever",
          value: "snork",
        },
      })
    ).toEqual([undefined]);
    expect(
      _whereCacheGet(whereCache, schema, "X", {
        where: {
          field: "var",
          column: "foo",
          value: "snork",
        },
      })
    ).toEqual([undefined]);
  });

  test("whereCacheCreate, where query", () => {
    const query = {
      where: {
        field: "queryField",
        column: "queryColumn",
        value: "queryValue",
      },
    };
    const wc = _whereCacheCreate(
      "field",
      {
        where: {
          field: "queryField",
          column: "queryColumn",
          value: "queryValue",
        },
      },
      [0, 1, 2]
    );
    expect(wc).toBeDefined();
    expect(wc).toEqual(
      expect.objectContaining({
        where: {
          field: {
            queryField: expect.any(Map),
          },
        },
      })
    );
    expect(wc.where.field.queryField.has("queryColumn")).toEqual(true);
    expect(wc.where.field.queryField.get("queryColumn")).toBeInstanceOf(Map);
    expect(
      wc.where.field.queryField.get("queryColumn").has("queryValue")
    ).toEqual(true);
    expect(_whereCacheGet(wc, schema, "field", query)).toEqual([0, 1, 2]);
  });

  test("whereCacheCreate, summarize query", () => {
    const query = {
      summarize: {
        method: "method",
        field: "queryField",
        column: "queryColumn",
        values: ["queryValue"],
      },
    };
    const wc = _whereCacheCreate("field", query, [0, 1, 2]);
    expect(_whereCacheGet(wc, schema, "field", query)).toEqual([0, 1, 2]);
  });

  test("whereCacheCreate, unknown query type", () => {
    expect(_whereCacheCreate("field", { foobar: true }, [1])).toEqual({});
  });

  test("whereCacheMerge, where queries", () => {
    let wc;

    // remember, will mutate dst
    const src = _whereCacheCreate(
      "field",
      { where: { field: "queryField", column: "queryColumn", value: "foo" } },
      ["foo"]
    );

    const dst1 = _whereCacheCreate(
      "field",
      { where: { field: "queryField", column: "queryColumn", value: "bar" } },
      ["dst1"]
    );
    wc = _whereCacheMerge(dst1, src);
    expect(
      _whereCacheGet(wc, schema, "field", {
        where: {
          field: "queryField",
          column: "queryColumn",
          value: "foo",
        },
      })
    ).toEqual(["foo"]);
    expect(
      _whereCacheGet(wc, schema, "field", {
        where: {
          field: "queryField",
          column: "queryColumn",
          value: "bar",
        },
      })
    ).toEqual(["dst1"]);

    const dst2 = _whereCacheCreate(
      "field",
      { where: { field: "queryField", column: "queryColumn", value: "bar" } },
      ["dst2"]
    );
    wc = _whereCacheMerge(dst2, dst1, src);
    expect(
      _whereCacheGet(wc, schema, "field", {
        where: {
          field: "queryField",
          column: "queryColumn",
          value: "foo",
        },
      })
    ).toEqual(["foo"]);
    expect(
      _whereCacheGet(wc, schema, "field", {
        where: {
          field: "queryField",
          column: "queryColumn",
          value: "bar",
        },
      })
    ).toEqual(["dst1"]);

    wc = _whereCacheMerge({}, src);
    expect(wc).toEqual(src);

    wc = _whereCacheMerge({ where: { field: { queryField: new Map() } } }, src);
    expect(wc).toEqual(src);
  });

  test("whereCacheMerge, mixed queries", () => {
    const wc = _whereCacheMerge(
      _whereCacheCreate(
        "field",
        {
          where: {
            field: "queryField",
            column: "queryColumn",
            value: "foo",
          },
        },
        ["a"]
      ),
      _whereCacheCreate(
        "field",
        {
          summarize: {
            method: "mean",
            field: "queryField",
            column: "queryColumn",
            values: ["foo", "bar", "baz"],
          },
        },
        ["b"]
      )
    );

    expect(
      _whereCacheGet(wc, schema, "field", {
        where: {
          field: "queryField",
          column: "queryColumn",
          value: "foo",
        },
      })
    ).toEqual(["a"]);

    expect(
      _whereCacheGet(wc, schema, "field", {
        summarize: {
          method: "mean",
          field: "queryField",
          column: "queryColumn",
          values: ["foo", "bar", "baz"],
        },
      })
    ).toEqual(["b"]);

    expect(
      _whereCacheGet(wc, schema, "field", {
        where: {
          field: "queryField",
          column: "queryColumn",
          value: "does-not-exist",
        },
      })
    ).toEqual([undefined]);

    expect(
      _whereCacheGet(wc, schema, "field", {
        summarize: {
          method: "no-such-method",
          field: "queryField",
          column: "queryColumn",
          values: ["does-not-exist"],
        },
      })
    ).toEqual([undefined]);
  });
});
