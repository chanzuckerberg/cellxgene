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
import { rangeFill } from "../../../src/util/range";

enableFetchMocks();

describe("AnnoMatrixCrossfilter", () => {
  let annoMatrix;
  let crossfilter;

  beforeEach(async () => {
    fetch.resetMocks(); // reset all fetch mocking state
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

      fetch.once(serverMocks.dataframeResponse(["louvain"], [obsLouvain]));
      let newCrossfilter = await crossfilter.select("obs", "louvain", {
        mode: "none",
      });

      expect(newCrossfilter.countSelected()).toEqual(0);
      expect(
        newCrossfilter.obsCrossfilter.hasDimension("obs/louvain")
      ).toBeTruthy();
      expect(fetch.mock.calls).toHaveLength(1);

      newCrossfilter = await crossfilter.select("obs", "louvain", {
        mode: "all",
      });
      expect(newCrossfilter.countSelected()).toEqual(annoMatrix.nObs);
    });

    test("simple column select", async () => {
      let xfltr;

      fetch.once(serverMocks.dataframeResponse(["louvain"], [obsLouvain]));
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
        (val, idx) => !["NK cells", "B cells"].includes(val) !== !selected[idx]
      );

      fetch.once(
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
      fetch.once(
        serverMocks.dataframeResponse(
          ["TEST"],
          [rangeFill(new Float32Array(nObs), 0, 0.1)]
        )
      );

      const xfltr = await crossfilter.select(
        "X",
        {
          field: "var",
          column: varIndex,
          value: "TYMP",
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
        field: "var",
        column: varIndex,
        value: "TYMP",
      });
      const values = df.icol(0).asArray();
      const selected = xfltr.allSelectedMask();
      values.every((val, idx) => !(val >= 0 && val <= 50) !== !selected[idx]);
      expect(selected.reduce((acc, val) => (val ? acc + 1 : acc), 0)).toEqual(
        xfltr.countSelected()
      );
    });

    test("spatial column select", async () => {
      fetch.once(
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
        mask[i] = true;
      }
      const annoMatrixSubset = isubsetMask(annoMatrix, mask);
      expect(annoMatrixSubset.nObs).toEqual(Math.floor(annoMatrix.nObs / 2));

      let xfltr = new AnnoMatrixObsCrossfilter(annoMatrixSubset);
      expect(xfltr.countSelected()).toEqual(annoMatrixSubset.nObs);

      fetch.once(serverMocks.dataframeResponse(["louvain"], [obsLouvain]));
      xfltr = await xfltr.select("obs", "louvain", {
        mode: "exact",
        values: ["NK cells", "B cells"],
      });

      expect(xfltr).toBeDefined();
      expect(xfltr.countSelected()).toEqual(240);

      const df = await annoMatrixSubset.fetch("obs", "louvain");
      const values = df.col("louvain").asArray();
      const selected = xfltr.allSelectedMask();
      values.every(
        (val, idx) => !["NK cells", "B cells"].includes(val) !== !selected[idx]
      );
    });

    test("select catches errors", async () => {
      await expect(crossfilter.select("NADA", "foo")).rejects.toThrow(
        "Unknown field name"
      );

      await expect(crossfilter.select("var", "foo")).rejects.toThrow(
        "unable to obsSelect upon the var dimension"
      );

      fetch.mockRejectOnce(new Error("unknown column name"));
      await expect(crossfilter.select("obs", "foo")).rejects.toThrow(
        "unknown column name"
      );
    });
  });

  describe("mutate matrix", () => {
    /*
    test the matrix mutators via crossfilter proxy
    */

    test("addObsColumn", async () => {
      expect(crossfilter.countSelected()).toBe(annoMatrix.nObs);
      expect(
        crossfilter.annoMatrix.getMatrixColumns("obs").includes("foo")
      ).toBeFalsy();
      const xfltr = crossfilter.addObsColumn(
        {
          name: "foo",
          type: "string",
        },
        Array,
        "toasty"
      );

      expect(xfltr.countSelected()).toBe(annoMatrix.nObs);
      expect(
        xfltr.annoMatrix.getMatrixColumns("obs").includes("foo")
      ).toBeTruthy();
      const df = await xfltr.annoMatrix.fetch("obs", "foo");
      expect(
        df
          .col("foo")
          .asArray()
          .every((v) => v === "toasty")
      ).toBeTruthy();
    });

    /*
    To be written - tests for:
      dropObsColumn
      renameObsColumn
      addObsAnnoCategory
      removeObsAnnoCategory
      setObsColumnValues
      resetObsColumnValues
    */
  });
});
