import { postUserErrorToast } from "../components/framework/toasters";
/*
Action creators for gene sets

Primarily used to keep the crossfilter and underlying data in sync with the UI.

The behavior manifest in these action creators:

    Delete a gene set, will
      * drop index & clear selection state on the gene set summary
      * drop index & clear selection state of each gene in the geneset

    Delete a gene from a gene set, will:
      * drop index & clear selection state on the gene set summary
      * drop index & clear selection state on the gene

    Add a gene to a gene set, will:
      * drop index & clear selection state on the gene set summary
      * will NOT touch the selection state for the gene

Note that crossfilter indices are lazy created, as needed.
*/

// eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any -- - FIXME: disabled temporarily on migrate to TS.
export const genesetDelete = (genesetName: any) => (
  // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any -- - FIXME: disabled temporarily on migrate to TS.
  dispatch: any,
  // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any -- - FIXME: disabled temporarily on migrate to TS.
  getState: any
) => {
  const state = getState();
  const { genesets } = state;
  const gs = genesets?.genesets?.get(genesetName) ?? {};
  const geneSymbols = Array.from(gs.genes.keys());
  const obsCrossfilter = dropGeneset(dispatch, state, genesetName, geneSymbols);
  if (genesetName === state.colors.colorAccessor) {
    dispatch({
      type: "reset colorscale",
    });
  }
  dispatch({
    type: "geneset: delete",
    genesetName,
    obsCrossfilter,
    annoMatrix: obsCrossfilter.annoMatrix,
  });
};

// eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any -- - FIXME: disabled temporarily on migrate to TS.
export const genesetAddGenes = (genesetName: any, genes: any) => async (
  // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any -- - FIXME: disabled temporarily on migrate to TS.
  dispatch: any,
  // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any -- - FIXME: disabled temporarily on migrate to TS.
  getState: any
) => {
  const state = getState();
  const { obsCrossfilter: prevObsCrossfilter, annoMatrix } = state;
  const { schema } = annoMatrix;
  const varIndex = schema.annotations.var.index;
  const df = await annoMatrix.fetch("var", varIndex);
  const geneNames = df.col(varIndex).asArray();
  // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
  genes = genes.reduce((acc: any, gene: any) => {
    if (geneNames.indexOf(gene.geneSymbol) === -1) {
      postUserErrorToast(
        `${gene.geneSymbol} doesn't appear to be a valid gene name.`
      );
    } else acc.push(gene);
    return acc;
  }, []);

  const obsCrossfilter = dropGenesetSummaryDimension(
    prevObsCrossfilter,
    state,
    genesetName
  );
  dispatch({
    type: "continuous metadata histogram cancel",
    continuousNamespace: { isGeneSetSummary: true },
    selection: genesetName,
  });
  return dispatch({
    type: "geneset: add genes",
    genesetName,
    genes,
    obsCrossfilter,
    annoMatrix: obsCrossfilter.annoMatrix,
  });
};

// eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any -- - FIXME: disabled temporarily on migrate to TS.
export const genesetDeleteGenes = (genesetName: any, geneSymbols: any) => (
  // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any -- - FIXME: disabled temporarily on migrate to TS.
  dispatch: any,
  // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any -- - FIXME: disabled temporarily on migrate to TS.
  getState: any
) => {
  const state = getState();
  const obsCrossfilter = dropGeneset(dispatch, state, genesetName, geneSymbols);
  return dispatch({
    type: "geneset: delete genes",
    genesetName,
    geneSymbols,
    obsCrossfilter,
    annoMatrix: obsCrossfilter.annoMatrix,
  });
};

/*
Private
*/

function dropGenesetSummaryDimension(
  // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
  obsCrossfilter: any,
  // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
  state: any,
  // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
  genesetName: any
) {
  const { annoMatrix, genesets } = state;
  const varIndex = annoMatrix.schema.annotations?.var?.index;
  const gs = genesets?.genesets?.get(genesetName) ?? {};
  const genes = Array.from(gs.genes.keys());
  const query = {
    summarize: {
      method: "mean",
      field: "var",
      column: varIndex,
      values: genes,
    },
  };
  return obsCrossfilter.dropDimension("X", query);
}

// eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
function dropGeneDimension(obsCrossfilter: any, state: any, gene: any) {
  const { annoMatrix } = state;
  const varIndex = annoMatrix.schema.annotations?.var?.index;
  const query = {
    where: {
      field: "var",
      column: varIndex,
      value: gene,
    },
  };
  return obsCrossfilter.dropDimension("X", query);
}

function dropGeneset(
  // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
  dispatch: any,
  // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
  state: any,
  // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
  genesetName: any,
  // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
  geneSymbols: any
) {
  const { obsCrossfilter: prevObsCrossfilter } = state;
  const obsCrossfilter = geneSymbols.reduce(
    // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
    (crossfilter: any, gene: any) =>
      dropGeneDimension(crossfilter, state, gene),
    dropGenesetSummaryDimension(prevObsCrossfilter, state, genesetName)
  );
  dispatch({
    type: "continuous metadata histogram cancel",
    continuousNamespace: { isGeneSetSummary: true },
    selection: genesetName,
  });
  // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
  geneSymbols.forEach((g: any) =>
    dispatch({
      type: "continuous metadata histogram cancel",
      continuousNamespace: { isUserDefined: true },
      selection: g,
    })
  );
  return obsCrossfilter;
}
