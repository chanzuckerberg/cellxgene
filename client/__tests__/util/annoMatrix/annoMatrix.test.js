import {
  AnnoMatrixLoader,
  AnnoMatrixObsCrossfilter,
} from "../../../src/util/annoMatrix/";
import { doJsonRequest } from "../../../src/util/actionHelpers";
import { Dataframe } from "../../../src/util/dataframe";

const baseDataUrl = "http://localhost:5005/api/v0.2";

const fetchJson = (path) => doJsonRequest(`${baseDataUrl}${path}`);
// const fetchBinary = (path) => doBinaryRequest(`${baseDataUrl}${path}`);

describe.skip("AnnoMatrix", () => {
  let annoMatrix;

  beforeEach(async () => {
    const schema = await fetchJson("/schema");
    annoMatrix = new AnnoMatrixLoader(baseDataUrl, schema.schema);
  });

  test("simple fetch", async () => {
    await expect(annoMatrix.fetch("obs", "name_0")).resolves.toBeInstanceOf(
      Dataframe
    );
    await expect(
      annoMatrix.fetch("obs", ["name_0", "n_genes"])
    ).resolves.toBeInstanceOf(Dataframe);
  });

  test("fetch all fields", async () => {
    ["emb", "var", "obs"].forEach(async (field) => {
      const columnNames = annoMatrix.getMatrixColumns(field);
      await expect(
        annoMatrix.fetch(field, columnNames.slice(2))
      ).resolves.toBeInstanceOf(Dataframe);
    });
  });

  test("fetch - test all query forms", async () => {
    // single string is a column name
    await expect(annoMatrix.fetch("obs", "n_genes")).resolves.toBeInstanceOf(
      Dataframe
    );

    // array of column names
    await expect(
      annoMatrix.fetch("obs", ["n_genes", "percent_mito"])
    ).resolves.toBeInstanceOf(Dataframe);

    // more complex value filter query, enumerated
    await expect(
      annoMatrix.fetch("X", {
        field: "var",
        column: annoMatrix.schema.annotations.var.index,
        value: "TYMP",
      })
    ).resolves.toBeInstanceOf(Dataframe);

    // more complex value filter query, range
    await expect(
      annoMatrix.fetch("X", [
        {
          field: "var",
          column: annoMatrix.schema.annotations.var.index,
          value: "SUMO3",
        },
        {
          field: "var",
          column: annoMatrix.schema.annotations.var.index,
          value: "TYMP",
        },
      ])
    ).resolves.toBeInstanceOf(Dataframe);
  });

  test("push and pop views", async () => {
    const am1 = annoMatrix.clip(0.1, 0.9);
    expect(am1.viewOf).toBe(annoMatrix);
    expect(am1.nObs).toEqual(annoMatrix.nObs);
    expect(am1.nVar).toEqual(annoMatrix.nVar);
    expect(am1.rowIndex).toBe(annoMatrix.rowIndex);

    const am2 = annoMatrix.clip(0.1, 0.9);
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

    const am1 = annoMatrix.isubset(rowList);
    const am2 = annoMatrix.isubsetMask(rowMask);
    expect(am1).not.toBe(am2);
    expect(am1.nObs).toEqual(2);
    expect(am1.nObs).toEqual(am2.nObs);
    expect(am1.nVar).toEqual(am2.nVar);

    const ng1 = await am1.fetch("obs", "n_genes");
    const ng2 = await am2.fetch("obs", "n_genes");
    expect(ng1.length).toEqual(ng2.length);
    expect(ng1.colIndex.labels()).toEqual(ng2.colIndex.labels());
    expect(ng1.col("n_genes").asArray()).toEqual(ng2.col("n_genes").asArray());
  });

  describe("add/drop column", () => {
    async function addDrop(base) {
      expect(base.getMatrixColumns("obs")).not.toContain("foo");
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
      const foo = await am1.fetch("obs", "foo");
      expect(foo).toBeDefined();
      expect(foo).toBeInstanceOf(Dataframe);
      expect(foo.length).toEqual(am1.nObs);
      expect(foo.col("foo").asArray()).toEqual(
        new Float32Array(am1.nObs).fill(0)
      );

      /* drop */
      const am2 = am1.dropObsColumn("foo");
      expect(base.getMatrixColumns("obs")).not.toContain("foo");
      expect(am2.getMatrixColumns("obs")).not.toContain("foo");
      await expect(am2.fetch("obs", "foo")).rejects.toThrow(
        "unknown column name"
      );
    }

    test("add/drop column, without view", async () => {
      await addDrop(annoMatrix);
    });

    test("add/drop column, with view", async () => {
      const am1 = annoMatrix.clip(0.1, 0.9);
      await addDrop(am1);

      const am2 = am1.isubset([0, 1, 2, 20, 30, 400]);
      await addDrop(am2);

      const am3 = annoMatrix.isubset([10, 0, 7, 3]);
      await addDrop(am3);

      const am4 = am3.clip(0, 1);
      await addDrop(am4);

      await am1.fetch("obs", am1.getMatrixColumns("obs"));
      await am2.fetch("obs", am2.getMatrixColumns("obs"));
      await am3.fetch("obs", am3.getMatrixColumns("obs"));
      await am4.fetch("obs", am4.getMatrixColumns("obs"));

      await addDrop(am1);
      await addDrop(am2);
      await addDrop(am3);
      await addDrop(am4);
    });
  });

  describe("setObsColumnValues", () => {
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
      const am1 = am.setObsColumnValues("test", whichRows, "yo");
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
      am = am1.dropObsColumn("test");
      await expect(am.fetch("obs", "test")).rejects.toThrow(
        "unknown column name"
      );
    }

    test("set, without a view", async () => {
      await addSetDrop(annoMatrix);
    });

    test("set, with a view", async () => {
      const am1 = annoMatrix.clip(0.1, 0.9);
      await addSetDrop(am1);

      const am2 = am1.isubset([0, 1, 2, 10, 20, 30, 400]);
      await addSetDrop(am2);

      const am3 = annoMatrix.isubset([10, 1, 0, 30, 2]);
      await addSetDrop(am3);

      await am1.fetch("obs", am1.getMatrixColumns("obs"));
      await am2.fetch("obs", am2.getMatrixColumns("obs"));
      await am3.fetch("obs", am3.getMatrixColumns("obs"));

      await addSetDrop(am1);
      await addSetDrop(am2);
      await addSetDrop(am3);
    });
  });
});

describe.skip("AnnoMatrixCrossfilter", () => {
  let annoMatrix;
  let crossfilter;

  beforeEach(async () => {
    const schema = await fetchJson("/schema");
    annoMatrix = new AnnoMatrixLoader(baseDataUrl, schema.schema);
    crossfilter = new AnnoMatrixObsCrossfilter(annoMatrix);
  });

  test("initial state of crossfilter", () => {
    expect(crossfilter).toBeDefined();
    expect(crossfilter.countSelectedObs()).toEqual(annoMatrix.nObs);
    expect(crossfilter.allSelectedMask()).toHaveLength(annoMatrix.nObs);
    expect(crossfilter.allSelectedMask()).toEqual(
      new Uint8Array(annoMatrix.nObs).fill(1)
    );
  });

  test("simple column select", async () => {
    let xfltr;
    xfltr = await crossfilter.selectObs("obs", "louvain", {
      mode: "exact",
      values: ["NK cells", "B cells"],
    });

    expect(xfltr).toBeDefined();
    expect(xfltr.countSelectedObs()).toEqual(496);

    const df = await annoMatrix.fetch("obs", "louvain");
    const values = df.col("louvain").asArray();
    const selected = xfltr.allSelectedMask();
    values.every(
      (val, idx) => !["NK cells", "B cells"].includes(val) !== !selected[idx]
    );

    xfltr = await xfltr.selectObs("obs", "n_genes", {
      mode: "range",
      lo: 0,
      hi: 500,
      inclusive: false,
    });
    expect(xfltr.countSelectedObs()).toEqual(33);
  });

  test("complex column select", async () => {
    let xfltr;
    xfltr = await crossfilter.selectObs(
      "X",
      {
        field: "var",
        column: annoMatrix.schema.annotations.var.index,
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
    expect(xfltr.countSelectedObs()).toEqual(887);

    const df = await annoMatrix.fetch("X", {
      field: "var",
      column: annoMatrix.schema.annotations.var.index,
      value: "TYMP",
    });
    const values = df.icol(0).asArray();
    const selected = xfltr.allSelectedMask();
    values.every((val, idx) => !(val >= 0 && val <= 50) !== !selected[idx]);
    expect(selected.reduce((acc, val) => (val ? acc + 1 : acc), 0)).toEqual(
      xfltr.countSelectedObs()
    );
  });

  test("select on subset", async () => {
    const mask = new Uint8Array(annoMatrix.nObs).fill(0);
    for (let i = 0; i < mask.length; i += 2) {
      mask[i] = true;
    }
    const annoMatrixSubset = annoMatrix.isubsetMask(mask);
    expect(annoMatrixSubset.nObs).toEqual(Math.floor(annoMatrix.nObs / 2));

    let xfltr = new AnnoMatrixObsCrossfilter(annoMatrixSubset);
    expect(xfltr.countSelectedObs()).toEqual(annoMatrixSubset.nObs);
    xfltr = await xfltr.selectObs("obs", "louvain", {
      mode: "exact",
      values: ["NK cells", "B cells"],
    });

    expect(xfltr).toBeDefined();
    expect(xfltr.countSelectedObs()).toEqual(240);

    const df = await annoMatrixSubset.fetch("obs", "louvain");
    const values = df.col("louvain").asArray();
    const selected = xfltr.allSelectedMask();
    values.every(
      (val, idx) => !["NK cells", "B cells"].includes(val) !== !selected[idx]
    );
  });
});
