import {
  _whereCacheGet,
  _whereCacheCreate,
  _whereCacheMerge,
} from "../../../src/annoMatrix/whereCache";

const schema = {};

describe("whereCache", () => {
  test("whereCacheGet - missing cache values", () => {
    expect(
      _whereCacheGet({}, schema, "X", {
        field: "var",
        column: "foo",
        value: "bar",
      })
    ).toEqual([undefined]);
    expect(
      _whereCacheGet({ X: {} }, schema, "X", {
        field: "var",
        column: "foo",
        value: "bar",
      })
    ).toEqual([undefined]);
    expect(
      _whereCacheGet({ X: { var: new Map() } }, schema, "X", {
        field: "var",
        column: "foo",
        value: "bar",
      })
    ).toEqual([undefined]);
    expect(
      _whereCacheGet(
        { X: { var: new Map([["foo", new Map()]]) } },
        schema,
        "X",
        {
          field: "var",
          column: "foo",
          value: "bar",
        }
      )
    ).toEqual([undefined]);
  });

  test("whereCacheGet - varied lookups", () => {
    const whereCache = {
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
    };

    expect(
      _whereCacheGet(whereCache, schema, "X", {
        field: "var",
        column: "foo",
        value: "bar",
      })
    ).toEqual([0]);
    expect(
      _whereCacheGet(whereCache, schema, "X", {
        field: "var",
        column: "foo",
        value: "baz",
      })
    ).toEqual([1, 2]);
    expect(_whereCacheGet(whereCache, schema, "Y", {})).toEqual([undefined]);
    expect(
      _whereCacheGet(whereCache, schema, "X", {
        field: "whoknows",
        column: "whatever",
        value: "snork",
      })
    ).toEqual([undefined]);
    expect(
      _whereCacheGet(whereCache, schema, "X", {
        field: "var",
        column: "whatever",
        value: "snork",
      })
    ).toEqual([undefined]);
    expect(
      _whereCacheGet(whereCache, schema, "X", {
        field: "var",
        column: "foo",
        value: "snork",
      })
    ).toEqual([undefined]);
  });

  test("whereCacheCreate", () => {
    const query = {
      field: "queryField",
      column: "queryColumn",
      value: "queryValue",
    };
    const wc = _whereCacheCreate(
      "field",
      { field: "queryField", column: "queryColumn", value: "queryValue" },
      [0, 1, 2]
    );
    expect(wc).toBeDefined();
    expect(wc).toEqual(
      expect.objectContaining({
        field: {
          queryField: expect.any(Map),
        },
      })
    );
    expect(wc.field.queryField.has("queryColumn")).toEqual(true);
    expect(wc.field.queryField.get("queryColumn")).toBeInstanceOf(Map);
    expect(wc.field.queryField.get("queryColumn").has("queryValue")).toEqual(
      true
    );
    expect(_whereCacheGet(wc, schema, "field", query)).toEqual([0, 1, 2]);
  });

  test("whereCacheMerge", () => {
    let wc;

    // remember, will mutate dst
    const src = _whereCacheCreate(
      "field",
      { field: "queryField", column: "queryColumn", value: "foo" },
      ["foo"]
    );

    const dst1 = _whereCacheCreate(
      "field",
      { field: "queryField", column: "queryColumn", value: "bar" },
      ["dst1"]
    );
    wc = _whereCacheMerge(dst1, src);
    expect(
      _whereCacheGet(wc, schema, "field", {
        field: "queryField",
        column: "queryColumn",
        value: "foo",
      })
    ).toEqual(["foo"]);
    expect(
      _whereCacheGet(wc, schema, "field", {
        field: "queryField",
        column: "queryColumn",
        value: "bar",
      })
    ).toEqual(["dst1"]);

    const dst2 = _whereCacheCreate(
      "field",
      { field: "queryField", column: "queryColumn", value: "bar" },
      ["dst2"]
    );
    wc = _whereCacheMerge(dst2, dst1, src);
    expect(
      _whereCacheGet(wc, schema, "field", {
        field: "queryField",
        column: "queryColumn",
        value: "foo",
      })
    ).toEqual(["foo"]);
    expect(
      _whereCacheGet(wc, schema, "field", {
        field: "queryField",
        column: "queryColumn",
        value: "bar",
      })
    ).toEqual(["dst1"]);

    wc = _whereCacheMerge({}, src);
    expect(wc).toEqual(src);

    wc = _whereCacheMerge({ field: { queryField: new Map() } }, src);
    expect(wc).toEqual(src);
  });
});
