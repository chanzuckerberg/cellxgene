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

export const genesetDelete = (genesetName) => (dispatch, getState) => {
  const state = getState();
  const { genesets } = state;
  const gs = genesets?.genesets?.get(genesetName) ?? {};
  const geneSymbols = Array.from(gs.genes.keys());
  const obsCrossfilter = dropGeneset(dispatch, state, genesetName, geneSymbols);
  dispatch({
    type: "geneset: delete",
    genesetName,
    obsCrossfilter,
    annoMatrix: obsCrossfilter.annoMatrix,
  });
};

export const genesetAddGenes = (genesetName, genes) => (dispatch, getState) => {
  const state = getState();
  const { obsCrossfilter: prevObsCrossfilter } = state;
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

export const genesetDeleteGenes = (genesetName, geneSymbols) => (
  dispatch,
  getState
) => {
  const state = getState();
  const obsCrossfilter = dropGeneset(dispatch, state, genesetName, geneSymbols);
  // const { obsCrossfilter: prevObsCrossfilter } = state;
  // const obsCrossfilter = geneSymbols.reduce(
  //   (crossfilter, gene) => dropGeneDimension(crossfilter, state, gene),
  //   dropGenesetSummaryDimension(prevObsCrossfilter, state, genesetName)
  // );
  // dispatch({
  //   type: "continuous metadata histogram cancel",
  //   continuousNamespace: { isGeneSetSummary: true },
  //   selection: genesetName,
  // });
  // geneSymbols.forEach((g) =>
  //   dispatch({
  //     type: "continuous metadata histogram cancel",
  //     continuousNamespace: { isUserDefined: true },
  //     selection: g,
  //   })
  // );
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

function dropGenesetSummaryDimension(obsCrossfilter, state, genesetName) {
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

function dropGeneDimension(obsCrossfilter, state, gene) {
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

function dropGeneset(dispatch, state, genesetName, geneSymbols) {
  const { obsCrossfilter: prevObsCrossfilter } = state;
  const obsCrossfilter = geneSymbols.reduce(
    (crossfilter, gene) => dropGeneDimension(crossfilter, state, gene),
    dropGenesetSummaryDimension(prevObsCrossfilter, state, genesetName)
  );
  dispatch({
    type: "continuous metadata histogram cancel",
    continuousNamespace: { isGeneSetSummary: true },
    selection: genesetName,
  });
  geneSymbols.forEach((g) =>
    dispatch({
      type: "continuous metadata histogram cancel",
      continuousNamespace: { isUserDefined: true },
      selection: g,
    })
  );
  return obsCrossfilter;
}
