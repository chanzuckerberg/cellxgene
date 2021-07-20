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

export const genesetDelete = (genesetName: any) => (
  dispatch: any,
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

export const genesetAddGenes = (genesetName: any, genes: any) => async (
  dispatch: any,
  getState: any
) => {
  const state = getState();
  const { obsCrossfilter: prevObsCrossfilter, annoMatrix } = state;
  const { schema } = annoMatrix;
  const varIndex = schema.annotations.var.index;
  const df = await annoMatrix.fetch("var", varIndex);
  const geneNames = df.col(varIndex).asArray();
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

export const genesetDeleteGenes = (genesetName: any, geneSymbols: any) => (
  dispatch: any,
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
  obsCrossfilter: any,
  state: any,
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
  dispatch: any,
  state: any,
  genesetName: any,
  geneSymbols: any
) {
  const { obsCrossfilter: prevObsCrossfilter } = state;
  const obsCrossfilter = geneSymbols.reduce(
    (crossfilter: any, gene: any) =>
      dropGeneDimension(crossfilter, state, gene),
    dropGenesetSummaryDimension(prevObsCrossfilter, state, genesetName)
  );
  dispatch({
    type: "continuous metadata histogram cancel",
    continuousNamespace: { isGeneSetSummary: true },
    selection: genesetName,
  });
  geneSymbols.forEach((g: any) =>
    dispatch({
      type: "continuous metadata histogram cancel",
      continuousNamespace: { isUserDefined: true },
      selection: g,
    })
  );
  return obsCrossfilter;
}
