import sha1 from "sha1";
import {
  _whereCacheGet,
  _whereCacheCreate,
  _whereCacheMerge,
} from "../../../src/annoMatrix/whereCache";
import { Field, Schema } from "../../../src/common/types/schema";
import { Query } from "../../../src/annoMatrix/query";

const schema = {} as Schema;

describe("whereCache", () => {
  test("whereCacheGet - where query, missing cache values", () => {
    expect(
      _whereCacheGet({}, schema, Field.X, {
        where: {
          field: Field.var,
          column: "foo",
          value: "bar",
        },
      })
    ).toEqual([undefined]);
    expect(
      _whereCacheGet({}, schema, Field.X, {
        summarize: {
          method: "mean",
          field: Field.var,
          column: "foo",
          values: ["bar"],
        },
      })
    ).toEqual([undefined]);
    expect(
      _whereCacheGet({ where: { X: {} } }, schema, Field.X, {
        where: {
          field: "var",
          column: "foo",
          value: "bar",
        },
      })
    ).toEqual([undefined]);
    expect(
      _whereCacheGet({ where: { X: { var: new Map() } } }, schema, Field.X, {
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
        Field.X,
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
      _whereCacheGet({}, schema, Field.X, {
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
        Field.X,
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
      _whereCacheGet(whereCache, schema, Field.X, {
        where: {
          field: "var",
          column: "foo",
          value: "bar",
        },
      })
    ).toEqual([0]);
    expect(
      _whereCacheGet(whereCache, schema, Field.X, {
        summarize: {
          method: "mean",
          field: "var",
          column: "foo",
          values: ["bar"],
        },
      })
    ).toEqual([0]);
    expect(
      _whereCacheGet(whereCache, schema, Field.X, {
        where: {
          field: "var",
          column: "foo",
          value: "baz",
        },
      })
    ).toEqual([1, 2]);
    expect(
      _whereCacheGet(whereCache, schema, Field.X, {
        summarize: {
          method: "mean",
          field: "var",
          column: "foo",
          values: ["baz"],
        },
      })
    ).toEqual([1, 2]);
    expect(
      // Force invalid field value Y
      _whereCacheGet(whereCache, schema, "Y" as Field, {} as Query)
    ).toEqual([undefined]);
    expect(
      _whereCacheGet(whereCache, schema, Field.X, {
        where: {
          field: "whoknows",
          column: "whatever",
          value: "snork",
        },
      })
    ).toEqual([undefined]);
    expect(
      _whereCacheGet(whereCache, schema, Field.X, {
        where: {
          field: "var",
          column: "whatever",
          value: "snork",
        },
      })
    ).toEqual([undefined]);
    expect(
      _whereCacheGet(whereCache, schema, Field.X, {
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
      Field.obs,
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
          [Field.obs]: {
            queryField: expect.any(Map),
          },
        },
      })
    );
    // @ts-expect-error ts-migrate(2531) FIXME: Object is possibly 'null'.
    // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
    expect((wc.where as any)[Field.obs].queryField.has("queryColumn")).toEqual(true);
    expect(
      // @ts-expect-error ts-migrate(2531) FIXME: Object is possibly 'null'.
      // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
      (wc.where as any).obs.queryField.get("queryColumn")
    ).toBeInstanceOf(Map);
    expect(
      // @ts-expect-error ts-migrate(2531) FIXME: Object is possibly 'null'.
      // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
      (wc.where as any).obs.queryField.get("queryColumn").has("queryValue")
    ).toEqual(true);
    // eslint-disable-next-line @typescript-eslint/no-non-null-assertion --- assert wc to be non-null
    expect(_whereCacheGet(wc!, schema, Field.obs, query)).toEqual([0, 1, 2]);
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
    const wc = _whereCacheCreate(Field.obs, query, [0, 1, 2]);
    // eslint-disable-next-line @typescript-eslint/no-non-null-assertion --- assert wc to be non-null
    expect(_whereCacheGet(wc!, schema, Field.obs, query)).toEqual([0, 1, 2]);
  });

  test("whereCacheCreate, unknown query type", () => {
    // @ts-expect-error --- force invalid query {foobar: true}
    expect(_whereCacheCreate(Field.obs, { foobar: true }, [1])).toEqual({});
  });

  test("whereCacheMerge, where queries", () => {
    let wc;

    // remember, will mutate dst
    const src = _whereCacheCreate(
      Field.obs,
      { where: { field: "queryField", column: "queryColumn", value: "foo" } },
      ["foo"]
    );

    const dst1 = _whereCacheCreate(
      Field.obs,
      { where: { field: "queryField", column: "queryColumn", value: "bar" } },
      ["dst1"]
    );
    // eslint-disable-next-line @typescript-eslint/no-non-null-assertion --- assert dst1 and src to be non-null
    wc = _whereCacheMerge(dst1!, src!);
    expect(
      _whereCacheGet(wc, schema, Field.obs, {
        where: {
          field: "queryField",
          column: "queryColumn",
          value: "foo",
        },
      })
    ).toEqual(["foo"]);
    expect(
      _whereCacheGet(wc, schema, Field.obs, {
        where: {
          field: "queryField",
          column: "queryColumn",
          value: "bar",
        },
      })
    ).toEqual(["dst1"]);

    const dst2 = _whereCacheCreate(
      Field.obs,
      { where: { field: "queryField", column: "queryColumn", value: "bar" } },
      ["dst2"]
    );
    // eslint-disable-next-line @typescript-eslint/no-non-null-assertion --- assert dst2m, dst1 and src to be non-null
    wc = _whereCacheMerge(dst2!, dst1!, src!);
    expect(
      _whereCacheGet(wc, schema, Field.obs, {
        where: {
          field: "queryField",
          column: "queryColumn",
          value: "foo",
        },
      })
    ).toEqual(["foo"]);
    expect(
      _whereCacheGet(wc, schema, Field.obs, {
        where: {
          field: "queryField",
          column: "queryColumn",
          value: "bar",
        },
      })
    ).toEqual(["dst1"]);

    // eslint-disable-next-line @typescript-eslint/no-non-null-assertion --- assert src to be non-null
    wc = _whereCacheMerge({}, src!);
    expect(wc).toEqual(src);

    // eslint-disable-next-line @typescript-eslint/no-non-null-assertion --- assert src to be non-null
    wc = _whereCacheMerge({ where: { obs: { queryField: new Map() } } }, src!);
    expect(wc).toEqual(src);
  });

  test("whereCacheMerge, mixed queries", () => {
    const wc = _whereCacheMerge(
      // eslint-disable-next-line @typescript-eslint/no-non-null-assertion --- assert wc to be non-null
      _whereCacheCreate(
        Field.obs,
        {
          where: {
            field: "queryField",
            column: "queryColumn",
            value: "foo",
          },
        },
        ["a"]
      )!,
      // eslint-disable-next-line @typescript-eslint/no-non-null-assertion --- assert wc to be non-null
      _whereCacheCreate(
        Field.obs,
        {
          summarize: {
            method: "mean",
            field: "queryField",
            column: "queryColumn",
            values: ["foo", "bar", "baz"],
          },
        },
        ["b"]
      )!
    );

    expect(
      _whereCacheGet(wc, schema, Field.obs, {
        where: {
          field: "queryField",
          column: "queryColumn",
          value: "foo",
        },
      })
    ).toEqual(["a"]);

    expect(
      _whereCacheGet(wc, schema, Field.obs, {
        summarize: {
          method: "mean",
          field: "queryField",
          column: "queryColumn",
          values: ["foo", "bar", "baz"],
        },
      })
    ).toEqual(["b"]);

    expect(
      _whereCacheGet(wc, schema, Field.obs, {
        where: {
          field: "queryField",
          column: "queryColumn",
          value: "does-not-exist",
        },
      })
    ).toEqual([undefined]);

    expect(
      _whereCacheGet(wc, schema, Field.obs, {
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
