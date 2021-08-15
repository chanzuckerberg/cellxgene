// these TWO statements MUST be first in the file, before any other imports
import { enableFetchMocks } from "jest-fetch-mock";
import * as serverMocks from "./serverMocks";
// OK, continue on!

import obsLouvain from "./louvain.json";
import obsNGenes from "./n_genes.json";
import embUmap from "./umap.json";

import {
  AnnoMatrixLoader,
  AnnoMatrixObsCrossfilter,
  isubsetMask,
} from "../../../src/annoMatrix";
import { Dataframe } from "../../../src/util/dataframe";
import { rangeFill } from "../../../src/util/range";

enableFetchMocks();

describe("AnnoMatrixCrossfilter", () => {
  // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
  let annoMatrix: any;
  // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
  let crossfilter: any;

  beforeEach(async () => {
    // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
    (fetch as any).resetMocks(); // reset all fetch mocking state
    // reset all fetch mocking state
    annoMatrix = new AnnoMatrixLoader(
      serverMocks.baseDataURL,
      serverMocks.schema.schema
    );
    crossfilter = new AnnoMatrixObsCrossfilter(annoMatrix);
  });

  test("initial state of crossfilter", () => {
    const { nObs } = annoMatrix;

    expect(crossfilter).toBeDefined();
    expect(crossfilter.size()).toEqual(nObs);
    expect(crossfilter.annoMatrix).toBe(annoMatrix);

    /* by default, everything should be selected, even if no data in cache */
    expect(crossfilter.countSelected()).toEqual(nObs);
    expect(crossfilter.allSelectedLabels()).toEqual(
      rangeFill(new Int32Array(nObs))
    );
    expect(crossfilter.allSelectedMask()).toEqual(new Uint8Array(nObs).fill(1));
    expect(crossfilter.fillByIsSelected(new Uint8Array(nObs), 2, 1)).toEqual(
      new Uint8Array(nObs).fill(2)
    );
  });

  describe("select", () => {
    /*
    test the selection state via crossfilter proxy
    */

    test("select loads index", async () => {
      /*
      Select should transparently load/create dimension index.

      Internal dimension names are field/col:col:col..., eg,

        obs:louvain
        emb:umap_0:umap_1

      */
      expect(crossfilter.obsCrossfilter.dimensionNames()).toEqual([]);
      expect(
        crossfilter.obsCrossfilter.hasDimension("obs/louvain")
      ).toBeFalsy();

      // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
      (fetch as any).once(
        serverMocks.dataframeResponse(["louvain"], [obsLouvain])
      );
      let newCrossfilter = await crossfilter.select("obs", "louvain", {
        mode: "none",
      });

      expect(newCrossfilter.countSelected()).toEqual(0);
      expect(
        newCrossfilter.obsCrossfilter.hasDimension("obs/louvain")
      ).toBeTruthy();
      // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
      expect((fetch as any).mock.calls).toHaveLength(1);

      newCrossfilter = await crossfilter.select("obs", "louvain", {
        mode: "all",
      });
      expect(newCrossfilter.countSelected()).toEqual(annoMatrix.nObs);
    });

    test("simple column select", async () => {
      let xfltr;

      // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
      (fetch as any).once(
        serverMocks.dataframeResponse(["louvain"], [obsLouvain])
      );
      xfltr = await crossfilter.select("obs", "louvain", {
        mode: "exact",
        values: ["NK cells", "B cells"],
      });

      expect(xfltr).toBeDefined();
      expect(xfltr.countSelected()).toEqual(496);
      expect(xfltr.allSelectedMask()).toEqual(
        Uint8Array.from(
          obsLouvain.map((val) =>
            val === "NK cells" || val === "B cells" ? 1 : 0
          )
        )
      );
      expect(xfltr.allSelectedLabels()).toEqual(
        Int32Array.from(
          obsLouvain.reduce((acc, val, idx) => {
            // @ts-expect-error ts-migrate(2345) FIXME: Argument of type 'number' is not assignable to par... Remove this comment to see the full error message
            if (val === "NK cells" || val === "B cells") acc.push(idx);
            return acc;
          }, [])
        )
      );
      expect(
        xfltr.fillByIsSelected(new Uint8Array(annoMatrix.nObs), 3, 1)
      ).toEqual(
        Uint8Array.from(
          obsLouvain.map((val) =>
            val === "NK cells" || val === "B cells" ? 3 : 1
          )
        )
      );

      const df = await annoMatrix.fetch("obs", "louvain");
      const values = df.col("louvain").asArray();
      const selected = xfltr.allSelectedMask();
      values.every(
        // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
        (val: any, idx: any) =>
          !["NK cells", "B cells"].includes(val) !== !selected[idx]
      );

      // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
      (fetch as any).once(
        serverMocks.dataframeResponse(["n_genes"], [new Int32Array(obsNGenes)])
      );
      xfltr = await xfltr.select("obs", "n_genes", {
        mode: "range",
        lo: 0,
        hi: 500,
        inclusive: false,
      });
      expect(xfltr.countSelected()).toEqual(33);
      expect(xfltr.allSelectedLabels()).toEqual(
        Int32Array.from(
          obsNGenes.reduce((acc, val, idx) => {
            const louvain = obsLouvain[idx];
            if (
              val >= 0 &&
              val < 500 &&
              (louvain === "NK cells" || louvain === "B cells")
            )
              // @ts-expect-error ts-migrate(2345) FIXME: Argument of type 'number' is not assignable to par... Remove this comment to see the full error message
              acc.push(idx);
            return acc;
          }, [])
        )
      );

      xfltr = await xfltr.selectAll();
      expect(xfltr.countSelected()).toEqual(annoMatrix.nObs);
    });

    test("join column select", async () => {
      const varIndex = annoMatrix.schema.annotations.var.index;

      const { nObs } = annoMatrix.schema.dataframe;
      // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
      (fetch as any).once(
        serverMocks.dataframeResponse(
          ["TEST"],
          [rangeFill(new Float32Array(nObs), 0, 0.1)]
        )
      );

      const xfltr = await crossfilter.select(
        "X",
        {
          where: {
            field: "var",
            column: varIndex,
            value: "TYMP",
          },
        },
        {
          mode: "range",
          lo: 0,
          hi: 50,
          inclusive: true,
        }
      );

      expect(xfltr).toBeDefined();
      expect(xfltr.countSelected()).toEqual(501);

      const df = await annoMatrix.fetch("X", {
        where: {
          field: "var",
          column: varIndex,
          value: "TYMP",
        },
      });
      const values = df.icol(0).asArray();
      const selected = xfltr.allSelectedMask();
      values.every(
        // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
        (val: any, idx: any) => !(val >= 0 && val <= 50) !== !selected[idx]
      );
      expect(
        // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
        selected.reduce((acc: any, val: any) => (val ? acc + 1 : acc), 0)
      ).toEqual(xfltr.countSelected());
    });

    test("spatial column select", async () => {
      // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
      (fetch as any).once(
        serverMocks.dataframeResponse(
          ["umap_0", "umap_1"],
          [Float32Array.from(embUmap[0]), Float32Array.from(embUmap[1])]
        )
      );
      const xfltr = await crossfilter.select("emb", "umap", {
        mode: "within-rect",
        minX: 0,
        minY: 0,
        maxX: 0.5,
        maxY: 0.5,
      });
      expect(xfltr.countSelected()).toEqual(16);
    });

    test("select on subset", async () => {
      const mask = new Uint8Array(annoMatrix.nObs).fill(0);
      for (let i = 0; i < mask.length; i += 2) {
        // @ts-expect-error ts-migrate(2322) FIXME: Type 'boolean' is not assignable to type 'number'.
        mask[i] = true;
      }
      const annoMatrixSubset = isubsetMask(annoMatrix, mask);
      expect(annoMatrixSubset.nObs).toEqual(Math.floor(annoMatrix.nObs / 2));

      let xfltr = new AnnoMatrixObsCrossfilter(annoMatrixSubset);
      expect(xfltr.countSelected()).toEqual(annoMatrixSubset.nObs);

      // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
      (fetch as any).once(
        serverMocks.dataframeResponse(["louvain"], [obsLouvain])
      );
      xfltr = await xfltr.select("obs", "louvain", {
        mode: "exact",
        values: ["NK cells", "B cells"],
      });

      expect(xfltr).toBeDefined();
      expect(xfltr.countSelected()).toEqual(240);

      const df = await annoMatrixSubset.fetch("obs", "louvain");
      expect(df).toBeDefined();
      const values = (df as Dataframe).col("louvain").asArray();
      const selected = xfltr.allSelectedMask();
      values.every(
        // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
        (val: any, idx: any) =>
          !["NK cells", "B cells"].includes(val) !== !selected[idx]
      );
    });

    test("select catches errors", async () => {
      await expect(crossfilter.select("NADA", "foo")).rejects.toThrow(
        "Unknown field name"
      );

      await expect(crossfilter.select("var", "foo")).rejects.toThrow(
        "unable to obsSelect upon the var dimension"
      );

      // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
      (fetch as any).mockRejectOnce(new Error("unknown column name"));
      await expect(crossfilter.select("obs", "foo")).rejects.toThrow(
        "unknown column name"
      );
    });
  });

  describe("mutate matrix", () => {
    /*
    test the matrix mutators via crossfilter proxy
    */
    // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
    async function helperAddTestCol(cf: any, colName: any, colSchema = null) {
      expect(
        cf.annoMatrix.getMatrixColumns("obs").includes(colName)
      ).toBeFalsy();

      if (colSchema === null) {
        // @ts-expect-error ts-migrate(2322) FIXME: Type '{ name: any; type: string; categories: strin... Remove this comment to see the full error message
        colSchema = {
          name: colName,
          type: "categorical",
          categories: ["toasty"],
        };
      }
      // @ts-expect-error ts-migrate(2531) FIXME: Object is possibly 'null'.
      colSchema.name = colName;
      // @ts-expect-error ts-migrate(2531) FIXME: Object is possibly 'null'.
      const initValue = colSchema.categories[0];
      const xfltr = cf.addObsColumn(colSchema, Array, initValue);
      expect(
        xfltr.annoMatrix.schema.annotations.obs.columns.filter(
          // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
          (v: any) => v.name === colName
        )
      ).toHaveLength(1);
      const df = await xfltr.annoMatrix.fetch("obs", colName);
      expect(df.hasCol(colName)).toBeTruthy();
      return xfltr;
    }

    test("addObsColumn", async () => {
      expect(crossfilter.countSelected()).toBe(annoMatrix.nObs);
      expect(
        crossfilter.annoMatrix.getMatrixColumns("obs").includes("foo")
      ).toBeFalsy();
      const xfltr = crossfilter.addObsColumn(
        { name: "foo", type: "categorical", categories: ["A"] },
        Array,
        "A"
      );

      // check schema updates correctly.
      expect(xfltr.countSelected()).toBe(annoMatrix.nObs);
      expect(
        xfltr.annoMatrix.getMatrixColumns("obs").includes("foo")
      ).toBeTruthy();
      expect(xfltr.annoMatrix.schema.annotations.obsByName.foo).toMatchObject({
        name: "foo",
        type: "categorical",
      });
      expect(
        xfltr.annoMatrix.schema.annotations.obs.columns.filter(
          // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
          (v: any) => v.name === "foo"
        )
      ).toHaveLength(1);

      // check data update.
      const df = await xfltr.annoMatrix.fetch("obs", "foo");
      expect(
        df
          .col("foo")
          .asArray() // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
          .every((v: any) => v === "A")
      ).toBeTruthy();

      // check that we catch dups
      expect(() =>
        xfltr.addObsColumn(
          { name: "foo", type: "categorical" },
          Array,
          "toasty"
        )
      ).toThrow("column already exists");
      expect(() =>
        xfltr.addObsColumn(
          { name: "louvain", type: "categorical" },
          Array,
          "toasty"
        )
      ).toThrow("column already exists");
    });

    test("dropObsColumn", async () => {
      let xfltr;

      /* check that we catch attempt to drop readonly dimension */
      expect(() => crossfilter.dropObsColumn("louvain")).toThrow(
        "Unknown or readonly obs column"
      );
      /* non-existent column */
      expect(() => crossfilter.dropObsColumn("does-not-exist")).toThrow(
        "Unknown or readonly obs column"
      );

      // add a column, then drop it.
      xfltr = await helperAddTestCol(crossfilter, "foo");
      xfltr = xfltr.dropObsColumn("foo");
      expect(
        xfltr.annoMatrix.schema.annotations.obs.columns.filter(
          // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
          (v: any) => v.name === "foo"
        )
      ).toHaveLength(0);
      expect(xfltr.annoMatrix.schema.annotations.obsByName.foo).toBeUndefined();
      // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
      (fetch as any).mockRejectOnce(new Error("unknown column name"));
      await expect(xfltr.annoMatrix.fetch("obs", "foo")).rejects.toThrow(
        "unknown column name"
      );

      // now same, but ensure we have built an index before doing the drop
      xfltr = await helperAddTestCol(crossfilter, "bar");
      xfltr = await xfltr.select("obs", "bar", {
        mode: "exact",
        values: "whatever",
      });
      xfltr = xfltr.dropObsColumn("bar");

      // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
      (fetch as any).mockRejectOnce(new Error("unknown column name"));
      await expect(xfltr.select("obs", "bar", { mode: "all" })).rejects.toThrow(
        "unknown column name"
      );
    });

    test("renameObsColumn", async () => {
      let xfltr;

      /* catch attempts to rename non-existent or readonly columns */
      expect(() =>
        crossfilter.renameObsColumn("does-not-exist", "foo")
      ).toThrow("Unknown or readonly obs column");
      expect(() => crossfilter.renameObsColumn("louvain", "foo")).toThrow(
        "Unknown or readonly obs column"
      );

      // add a column, then rename it.
      xfltr = await helperAddTestCol(crossfilter, "foo");
      xfltr = xfltr.renameObsColumn("foo", "bar");
      expect(xfltr.annoMatrix.getColumnSchema("obs", "foo")).toBeUndefined();
      expect(xfltr.annoMatrix.getColumnSchema("obs", "bar")).toMatchObject({
        name: "bar",
        type: "categorical",
      });

      // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
      (fetch as any).mockRejectOnce(new Error("unknown column name"));
      await expect(xfltr.annoMatrix.fetch("obs", "foo")).rejects.toThrow(
        "unknown column name"
      );
      const df = await xfltr.annoMatrix.fetch("obs", "bar");
      expect(df.hasCol("bar")).toBeTruthy();

      // now same, but ensure we have built an index before doing the rename
      xfltr = await helperAddTestCol(crossfilter, "bar");
      xfltr = await xfltr.select("obs", "bar", {
        mode: "exact",
        values: "whatever",
      });
      xfltr = xfltr.renameObsColumn("bar", "xyz");

      // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
      (fetch as any).mockRejectOnce(new Error("unknown column name"));
      await expect(xfltr.select("obs", "bar", { mode: "all" })).rejects.toThrow(
        "unknown column name"
      );
      await expect(
        xfltr.select("obs", "xyz", { mode: "none" })
      ).resolves.toBeInstanceOf(AnnoMatrixObsCrossfilter);
    });

    test("addObsAnnoCategory", async () => {
      // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
      let xfltr: any;

      // catch unknown or readonly columns
      expect(() => crossfilter.addObsAnnoCategory("louvain", "mumble")).toThrow(
        "Unknown or readonly obs column"
      );
      expect(() =>
        crossfilter.addObsAnnoCategory("undefined-name", "mumble")
      ).toThrow("Unknown or readonly obs column");

      // add a column and then add category to it
      // @ts-expect-error ts-migrate(2345) FIXME: Argument of type '{ name: string; type: string; ca... Remove this comment to see the full error message
      xfltr = await helperAddTestCol(crossfilter, "foo", {
        name: "foo",
        type: "categorical",
        categories: ["unassigned"],
      });
      xfltr = xfltr.addObsAnnoCategory("foo", "a-new-label");
      expect(xfltr.annoMatrix.getColumnSchema("obs", "foo")).toMatchObject({
        name: "foo",
        type: "categorical",
        categories: expect.arrayContaining(["a-new-label", "unassigned"]),
      });

      // do it again, dup; should throw
      expect(() => xfltr.addObsAnnoCategory("foo", "a-new-label")).toThrow(
        "category already exists"
      );

      // now same, but ensure we have built an index before doing the operation
      // @ts-expect-error ts-migrate(2345) FIXME: Argument of type '{ name: string; type: string; ca... Remove this comment to see the full error message
      xfltr = await helperAddTestCol(crossfilter, "bar", {
        name: "bar",
        type: "categorical",
        categories: ["unassigned"],
      });
      xfltr = await xfltr.select("obs", "bar", {
        mode: "exact",
        values: "something",
      });
      xfltr = xfltr.addObsAnnoCategory("bar", "a-new-label");
      expect(xfltr.annoMatrix.getColumnSchema("obs", "bar")).toMatchObject({
        name: "bar",
        type: "categorical",
        categories: expect.arrayContaining(["a-new-label", "unassigned"]),
      });
    });

    test("removeObsAnnoCategory", async () => {
      let xfltr;

      // catch unknown or readonly categories
      await expect(() =>
        crossfilter.removeObsAnnoCategory("louvain", "mumble", "unassigned")
      ).rejects.toThrow("Unknown or readonly obs column");
      await expect(() =>
        crossfilter.removeObsAnnoCategory("undefined-name", "mumble")
      ).rejects.toThrow("Unknown or readonly obs column");

      // @ts-expect-error ts-migrate(2345) FIXME: Argument of type '{ name: string; type: string; ca... Remove this comment to see the full error message
      xfltr = await helperAddTestCol(crossfilter, "foo", {
        name: "foo",
        type: "categorical",
        categories: ["unassigned", "red", "green", "blue"],
      });
      xfltr = await xfltr.select("obs", "foo", { mode: "all" });
      expect(
        (await xfltr.annoMatrix.fetch("obs", "foo"))
          .col("foo")
          .asArray() // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
          .every((v: any) => v === "unassigned")
      ).toBeTruthy();
      expect(xfltr.annoMatrix.getColumnSchema("obs", "foo")).toMatchObject({
        name: "foo",
        type: "categorical",
        categories: expect.arrayContaining([
          "unassigned",
          "red",
          "green",
          "blue",
        ]),
      });

      // remove an unused category
      const xfltr1 = await xfltr.removeObsAnnoCategory("foo", "red", "mumble");
      expect(
        (await xfltr1.annoMatrix.fetch("obs", "foo"))
          .col("foo")
          .asArray() // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
          .every((v: any) => v === "unassigned")
      ).toBeTruthy();
      expect(xfltr1.annoMatrix.getColumnSchema("obs", "foo")).toMatchObject({
        name: "foo",
        type: "categorical",
        categories: expect.arrayContaining([
          "unassigned",
          "green",
          "blue",
          "mumble",
        ]),
      });

      // remove a used category
      const xfltr2 = await xfltr.removeObsAnnoCategory(
        "foo",
        "unassigned",
        "red"
      );
      expect(
        (await xfltr2.annoMatrix.fetch("obs", "foo"))
          .col("foo")
          .asArray() // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
          .every((v: any) => v === "red")
      ).toBeTruthy();
      expect(xfltr2.annoMatrix.getColumnSchema("obs", "foo")).toMatchObject({
        name: "foo",
        type: "categorical",
        categories: expect.arrayContaining(["green", "blue", "red"]),
      });
    });

    test("setObsColumnValues", async () => {
      // catch unknown or readonly categories
      await expect(() =>
        crossfilter.setObsColumnValues("louvain", [0, 1], "unassigned")
      ).rejects.toThrow("Unknown or readonly obs column");
      await expect(() =>
        crossfilter.setObsColumnValues("undefined-name", [0], "mumble")
      ).rejects.toThrow("Unknown or readonly obs column");

      // @ts-expect-error ts-migrate(2345) FIXME: Argument of type '{ name: string; type: string; ca... Remove this comment to see the full error message
      let xfltr = await helperAddTestCol(crossfilter, "foo", {
        name: "foo",
        type: "categorical",
        categories: ["unassigned", "red", "green", "blue"],
      });
      xfltr = await xfltr.select("obs", "foo", { mode: "all" });

      // catch unknown row label
      await expect(() =>
        xfltr.setObsColumnValues("foo", [-1], "red")
      ).rejects.toThrow("Unknown row label");

      // set a few rows
      expect(
        (await xfltr.annoMatrix.fetch("obs", "foo"))
          .col("foo")
          .asArray() // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
          .every((v: any) => v === "unassigned")
      ).toBeTruthy();
      const xfltr1 = await xfltr.setObsColumnValues("foo", [0, 10], "purple");
      expect(
        (await xfltr1.annoMatrix.fetch("obs", "foo"))
          .col("foo")
          .asArray()
          .every(
            // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
            (v: any, i: any) =>
              v === "unassigned" || (v === "purple" && (i === 0 || i === 10))
          )
      ).toBeTruthy();
      expect(xfltr1.annoMatrix.getColumnSchema("obs", "foo")).toMatchObject({
        name: "foo",
        type: "categorical",
        categories: expect.arrayContaining([
          "unassigned",
          "red",
          "green",
          "blue",
          "purple",
        ]),
      });

      expect(xfltr1.countSelected()).toEqual(xfltr1.annoMatrix.nObs);
      const xfltr2 = await xfltr1.select("obs", "foo", {
        mode: "exact",
        values: ["purple"],
      });
      expect(xfltr2.countSelected()).toEqual(2);
      expect(xfltr2.allSelectedLabels()).toEqual(Array.from([0, 10]));
    });

    test("resetObsColumnValues", async () => {
      // catch unknown or readonly categories
      await expect(() =>
        crossfilter.resetObsColumnValues("louvain", "red", "blue")
      ).rejects.toThrow("Unknown or readonly obs column");
      await expect(() =>
        crossfilter.resetObsColumnValues("undefined-name", "red", "blue")
      ).rejects.toThrow("Unknown or readonly obs column");

      // @ts-expect-error ts-migrate(2345) FIXME: Argument of type '{ name: string; type: string; ca... Remove this comment to see the full error message
      let xfltr = await helperAddTestCol(crossfilter, "foo", {
        name: "foo",
        type: "categorical",
        categories: ["unassigned", "red", "green", "blue"],
      });
      xfltr = await xfltr.select("obs", "foo", {
        mode: "exact",
        values: "red",
      });

      // catch unknown category name label
      await expect(() =>
        xfltr.resetObsColumnValues("foo", "unknown-label", "red")
      ).rejects.toThrow("unknown category");

      let xfltr1 = await xfltr.setObsColumnValues("foo", [0, 10], "purple");
      xfltr1 = await xfltr1.select("obs", "foo", {
        mode: "exact",
        values: "purple",
      });
      expect(
        (await xfltr1.annoMatrix.fetch("obs", "foo"))
          .col("foo")
          .asArray() // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
          .filter((v: any) => v === "purple")
      ).toHaveLength(2);

      xfltr1 = await xfltr1.resetObsColumnValues("foo", "purple", "magenta");
      expect(
        (await xfltr1.annoMatrix.fetch("obs", "foo"))
          .col("foo")
          .asArray() // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
          .filter((v: any) => v === "magenta")
      ).toHaveLength(2);
      expect(
        (await xfltr1.annoMatrix.fetch("obs", "foo"))
          .col("foo")
          .asArray() // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
          .filter((v: any) => v === "purple")
      ).toHaveLength(0);
      expect(xfltr1.annoMatrix.getColumnSchema("obs", "foo")).toMatchObject({
        name: "foo",
        type: "categorical",
        categories: expect.arrayContaining([
          "unassigned",
          "red",
          "green",
          "blue",
          "purple",
          "magenta",
        ]),
      });
    });
  });

  describe("edge cases", () => {
    test("transition from empty annoMatrix", async () => {
      // select before fetch needs to work
      // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
      (fetch as any).once(
        serverMocks.dataframeResponse(["louvain"], [obsLouvain])
      );
      const xfltr = await crossfilter.select("obs", "louvain", {
        mode: "exact",
        values: "B cells",
      });
      // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
      expect((fetch as any).mock.calls).toHaveLength(1);
      expect(xfltr.obsCrossfilter.hasDimension("obs/louvain")).toBeTruthy();
      expect(xfltr.obsCrossfilter.all()).toBe(xfltr.annoMatrix._cache.obs);
      expect(xfltr.countSelected()).toEqual(
        obsLouvain.reduce(
          (count, v) => (v === "B cells" ? count + 1 : count),
          0
        )
      );
    });
  });
});
