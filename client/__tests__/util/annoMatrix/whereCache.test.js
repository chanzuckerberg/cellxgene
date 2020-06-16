import {
  whereCacheGet,
  whereCacheCreate,
  whereCacheMerge,
} from "../../../src/util/annoMatrix/whereCache";

const schema = {};

describe("whereCache", () => {
  test("whereCacheGet - missing cache values", () => {
    expect(
      whereCacheGet({}, schema, "X", {
        field: "var",
        column: "foo",
        value: "bar",
      })
    ).toEqual([undefined]);
    expect(
      whereCacheGet({ X: {} }, schema, "X", {
        field: "var",
        column: "foo",
        value: "bar",
      })
    ).toEqual([undefined]);
    expect(
      whereCacheGet({ X: { var: new Map() } }, schema, "X", {
        field: "var",
        column: "foo",
        value: "bar",
      })
    ).toEqual([undefined]);
    expect(
      whereCacheGet(
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
      whereCacheGet(whereCache, schema, "X", {
        field: "var",
        column: "foo",
        value: "bar",
      })
    ).toEqual([0]);
    expect(
      whereCacheGet(whereCache, schema, "X", {
        field: "var",
        column: "foo",
        value: "baz",
      })
    ).toEqual([1, 2]);
    expect(whereCacheGet(whereCache, schema, "Y", {})).toEqual([undefined]);
    expect(
      whereCacheGet(whereCache, schema, "X", {
        field: "whoknows",
        column: "whatever",
        value: "snork",
      })
    ).toEqual([undefined]);
    expect(
      whereCacheGet(whereCache, schema, "X", {
        field: "var",
        column: "whatever",
        value: "snork",
      })
    ).toEqual([undefined]);
    expect(
      whereCacheGet(whereCache, schema, "X", {
        field: "var",
        column: "foo",
        value: "snork",
      })
    ).toEqual([undefined]);
  });

  test("whereCacheCreate", () => {
    let wc;
    let query;

    query = { field: "queryField", column: "queryColumn", value: "queryValue" };
    wc = whereCacheCreate(
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
    expect(whereCacheGet(wc, schema, "field", query)).toEqual([0, 1, 2]);
  });

  test("whereCacheMerge", () => {
    let wc;

    // remember, will mutate dst
    const src = whereCacheCreate(
      "field",
      { field: "queryField", column: "queryColumn", value: "foo" },
      ["foo"]
    );

    const dst1 = whereCacheCreate(
      "field",
      { field: "queryField", column: "queryColumn", value: "bar" },
      ["dst1"]
    );
    wc = whereCacheMerge(dst1, src);
    expect(
      whereCacheGet(wc, schema, "field", {
        field: "queryField",
        column: "queryColumn",
        value: "foo",
      })
    ).toEqual(["foo"]);
    expect(
      whereCacheGet(wc, schema, "field", {
        field: "queryField",
        column: "queryColumn",
        value: "bar",
      })
    ).toEqual(["dst1"]);

    const dst2 = whereCacheCreate(
      "field",
      { field: "queryField", column: "queryColumn", value: "bar" },
      ["dst2"]
    );
    wc = whereCacheMerge(dst2, dst1, src);
    expect(
      whereCacheGet(wc, schema, "field", {
        field: "queryField",
        column: "queryColumn",
        value: "foo",
      })
    ).toEqual(["foo"]);
    expect(
      whereCacheGet(wc, schema, "field", {
        field: "queryField",
        column: "queryColumn",
        value: "bar",
      })
    ).toEqual(["dst1"]);

    wc = whereCacheMerge({}, src);
    expect(wc).toEqual(src);

    wc = whereCacheMerge({ field: { queryField: new Map() } }, src);
    expect(wc).toEqual(src);
  });
});
