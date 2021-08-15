// these TWO statements MUST be first in the file, before any other imports
import { enableFetchMocks } from "jest-fetch-mock";
import * as serverMocks from "./serverMocks";
// OK, continue on!

import {
  AnnoMatrixLoader,
  clip,
  isubset,
  isubsetMask,
} from "../../../src/annoMatrix";
import { Dataframe } from "../../../src/util/dataframe";

enableFetchMocks();

describe("AnnoMatrix", () => {
  // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
  let annoMatrix: any;

  beforeEach(async () => {
    // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
    (fetch as any).resetMocks(); // reset all fetch mocking state
    // reset all fetch mocking state
    annoMatrix = new AnnoMatrixLoader(
      serverMocks.baseDataURL,
      serverMocks.schema.schema
    );
  });

  describe("basics", () => {
    test("annomatrix static checks", () => {
      expect(annoMatrix).toBeDefined();
      expect(annoMatrix.schema).toMatchObject(serverMocks.schema.schema);
      expect(annoMatrix.nObs).toEqual(serverMocks.schema.schema.dataframe.nObs);
      expect(annoMatrix.nVar).toEqual(serverMocks.schema.schema.dataframe.nVar);
      expect(annoMatrix.isView).toBeFalsy();
      expect(annoMatrix.viewOf).toBeUndefined();
      expect(annoMatrix.rowIndex).toBeDefined();
    });

    test("simple single column fetch", async () => {
      // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
      (fetch as any).once(serverMocks.annotationsObs(["name_0"]));

      const df = await annoMatrix.fetch("obs", "name_0");
      expect(df).toBeInstanceOf(Dataframe);
      expect(df.colIndex.labels()).toEqual(["name_0"]);
      expect(df.dims).toEqual([annoMatrix.nObs, 1]);
    });

    test("simple multi column fetch", async () => {
      // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
      (fetch as any)
        .once(serverMocks.annotationsObs(["name_0"]))
        .once(serverMocks.annotationsObs(["n_genes"]));

      await expect(
        annoMatrix.fetch("obs", ["name_0", "n_genes"])
      ).resolves.toBeInstanceOf(Dataframe);
    });

    describe("fetch from field", () => {
      // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
      const getLastTwo = async (field: any) => {
        const columnNames = annoMatrix.getMatrixColumns(field).slice(-2);
        // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
        (fetch as any).mockResponses(
          ...columnNames.map(() => serverMocks.responder)
        );
        await expect(
          annoMatrix.fetch(field, columnNames)
        ).resolves.toBeInstanceOf(Dataframe);
      };

      test("obs", async () => getLastTwo("obs"));
      test("var", async () => getLastTwo("var"));
      test("emb", async () => getLastTwo("emb"));
    });

    test("fetch - test all query forms", async () => {
      // single string is a column name
      // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
      (fetch as any).once(serverMocks.annotationsObs(["n_genes"]));
      await expect(annoMatrix.fetch("obs", "n_genes")).resolves.toBeInstanceOf(
        Dataframe
      );

      // array of column names, expecting n_genes to be cached.
      // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
      (fetch as any).once(serverMocks.annotationsObs(["percent_mito"]));
      await expect(
        annoMatrix.fetch("obs", ["n_genes", "percent_mito"])
      ).resolves.toBeInstanceOf(Dataframe);

      // more complex value filter query, enumerated
      // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
      (fetch as any).once(serverMocks.responder);
      await expect(
        annoMatrix.fetch("X", {
          where: {
            field: "var",
            column: annoMatrix.schema.annotations.var.index,
            value: "TYMP",
          },
        })
      ).resolves.toBeInstanceOf(Dataframe);

      // more complex value filter query, range
      const varIndex = annoMatrix.schema.annotations.var.index;
      // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
      (fetch as any)
        .once(
          serverMocks.withExpected("/data/var", [[`var:${varIndex}`, "SUMO3"]])
        )
        .once(
          serverMocks.withExpected("/data/var", [[`var:${varIndex}`, "TYMP"]])
        );
      await expect(
        annoMatrix.fetch("X", [
          {
            where: {
              field: "var",
              column: varIndex,
              value: "SUMO3",
            },
          },
          {
            where: {
              field: "var",
              column: varIndex,
              value: "TYMP",
            },
          },
        ])
      ).resolves.toBeInstanceOf(Dataframe);
      // XXX inspect the wherecache?
    });

    test("push and pop views", async () => {
      const am1 = clip(annoMatrix, 0.1, 0.9);
      expect(am1.viewOf).toBe(annoMatrix);
      expect(am1.nObs).toEqual(annoMatrix.nObs);
      expect(am1.nVar).toEqual(annoMatrix.nVar);
      expect(am1.rowIndex).toBe(annoMatrix.rowIndex);

      const am2 = clip(annoMatrix, 0.1, 0.9);
      expect(am2.viewOf).toBe(annoMatrix);
      expect(am2).not.toBe(am1);
      expect(am2.rowIndex).toBe(annoMatrix.rowIndex);
    });

    test("schema accessors", () => {
      expect(annoMatrix.getMatrixFields()).toEqual(
        expect.arrayContaining(["X", "obs", "emb", "var"])
      );
      expect(annoMatrix.getMatrixColumns("obs")).toEqual(
        expect.arrayContaining(["name_0", "n_genes", "louvain"])
      );
      expect(annoMatrix.getColumnSchema("emb", "umap")).toEqual({
        name: "umap",
        dims: ["umap_0", "umap_1"],
        type: "float32",
      });
      expect(annoMatrix.getColumnDimensions("emb", "umap")).toEqual([
        "umap_0",
        "umap_1",
      ]);
    });

    /*
  	test the mask & label access to subset via isubset and isubsetMask
  	*/
    test("isubset", async () => {
      const rowList = [0, 10];
      const rowMask = new Uint8Array(annoMatrix.nObs);
      for (let i = 0; i < rowList.length; i += 1) {
        rowMask[rowList[i]] = 1;
      }

      const am1 = isubset(annoMatrix, rowList);
      const am2 = isubsetMask(annoMatrix, rowMask);
      expect(am1).not.toBe(am2);
      expect(am1.nObs).toEqual(2);
      expect(am1.nObs).toEqual(am2.nObs);
      expect(am1.nVar).toEqual(am2.nVar);

      // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
      (fetch as any)
        .once(serverMocks.annotationsObs(["n_genes"]))
        .once(serverMocks.annotationsObs(["n_genes"]));
      const ng1 = (await am1.fetch("obs", "n_genes")) as Dataframe;
      const ng2 = (await am2.fetch("obs", "n_genes")) as Dataframe;
      expect(ng1).toBeDefined();
      expect(ng2).toBeDefined();
      expect(ng1).toHaveLength(ng2.length);
      expect(ng1.colIndex.labels()).toEqual(ng2.colIndex.labels());
      expect(ng1.col("n_genes").asArray()).toEqual(
        ng2.col("n_genes").asArray()
      );
    });
  });

  describe("add/drop column", () => {
    // @ts-expect-error ts-migrate(7006) FIXME: Parameter 'base' implicitly has an 'any' type.
    async function addDrop(base) {
      expect(base.getMatrixColumns("obs")).not.toContain("foo");
      // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
      (fetch as any).mockRejectOnce(new Error("unknown column name"));
      await expect(base.fetch("obs", "foo")).rejects.toThrow(
        "unknown column name"
      );

      /* add */
      const am1 = base.addObsColumn(
        { name: "foo", type: "float32", writable: true },
        Float32Array,
        0
      );
      expect(base.getMatrixColumns("obs")).not.toContain("foo");
      expect(am1.getMatrixColumns("obs")).toContain("foo");
      const foo: Dataframe = await am1.fetch("obs", "foo");
      expect(foo).toBeDefined();
      expect(foo).toBeInstanceOf(Dataframe);
      expect(foo).toHaveLength(am1.nObs);
      expect(foo.col("foo").asArray()).toEqual(
        new Float32Array(am1.nObs).fill(0)
      );

      /* drop */
      const am2 = am1.dropObsColumn("foo");
      expect(base.getMatrixColumns("obs")).not.toContain("foo");
      expect(am2.getMatrixColumns("obs")).not.toContain("foo");
      // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
      (fetch as any).mockRejectOnce(new Error("unknown column name"));
      await expect(am2.fetch("obs", "foo")).rejects.toThrow(
        "unknown column name"
      );
    }

    test("add/drop column, without view", async () => {
      await addDrop(annoMatrix);
    });

    test("add/drop column, with view", async () => {
      const am1 = clip(annoMatrix, 0.1, 0.9);
      await addDrop(am1);

      const am2 = isubset(am1, [0, 1, 2, 20, 30, 400]);
      await addDrop(am2);

      const am3 = isubset(annoMatrix, [10, 0, 7, 3]);
      await addDrop(am3);

      const am4 = clip(am3, 0, 1);
      await addDrop(am4);

      // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
      (fetch as any).mockResponse(serverMocks.responder);

      await am1.fetch("obs", am1.getMatrixColumns("obs"));
      await am2.fetch("obs", am2.getMatrixColumns("obs"));
      await am3.fetch("obs", am3.getMatrixColumns("obs"));
      await am4.fetch("obs", am4.getMatrixColumns("obs"));

      // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
      (fetch as any).resetMocks();

      await addDrop(am1);
      await addDrop(am2);
      await addDrop(am3);
      await addDrop(am4);
    });
  });

  describe("setObsColumnValues", () => {
    // @ts-expect-error ts-migrate(7006) FIXME: Parameter 'base' implicitly has an 'any' type.
    async function addSetDrop(base) {
      /* add column */
      let am = base.addObsColumn(
        {
          name: "test",
          type: "categorical",
          categories: ["unassigned", "red", "green"],
          writable: true,
        },
        Array,
        "unassigned"
      );

      const testVal = await am.fetch("obs", "test");
      expect(testVal.col("test").asArray()).toEqual(
        new Array(am.nObs).fill("unassigned")
      );

      /* set values in column */
      const whichRows = [1, 2, 10];
      const am1 = await am.setObsColumnValues("test", whichRows, "yo");
      const testVal1 = await am1.fetch("obs", "test");
      const expt = new Array(am1.nObs).fill("unassigned");
      for (let i = 0; i < whichRows.length; i += 1) {
        const offset = am1.rowIndex.getOffset(whichRows[i]);
        expt[offset] = "yo";
      }
      expect(testVal1).not.toBe(testVal);
      expect(testVal1.col("test").asArray()).toEqual(expt);
      expect(am1.getColumnSchema("obs", "test").type).toBe("categorical");
      expect(am1.getColumnSchema("obs", "test").categories).toEqual(
        expect.arrayContaining(["unassigned", "red", "green", "yo"])
      );

      /* drop column */
      // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
      (fetch as any).mockRejectOnce(new Error("unknown column name"));
      am = am1.dropObsColumn("test");
      await expect(am.fetch("obs", "test")).rejects.toThrow(
        "unknown column name"
      );
    }

    test("set, without a view", async () => {
      await addSetDrop(annoMatrix);
    });

    test("set, with a view", async () => {
      const am1 = clip(annoMatrix, 0.1, 0.9);
      await addSetDrop(am1);

      const am2 = isubset(am1, [0, 1, 2, 10, 20, 30, 400]);
      await addSetDrop(am2);

      const am3 = isubset(annoMatrix, [10, 1, 0, 30, 2]);
      await addSetDrop(am3);

      // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
      (fetch as any).mockResponse(serverMocks.responder);

      await am1.fetch("obs", am1.getMatrixColumns("obs"));
      await am2.fetch("obs", am2.getMatrixColumns("obs"));
      await am3.fetch("obs", am3.getMatrixColumns("obs"));

      await addSetDrop(am1);
      await addSetDrop(am2);
      await addSetDrop(am3);
    });
  });
});
