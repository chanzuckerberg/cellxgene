import { postUserErrorToast } from "../components/framework/toasters";
/*
Action creators for sample sets

Primarily used to keep the crossfilter and underlying data in sync with the UI.

The behavior manifest in these action creators:

    Delete a sample set, will
      * drop index & clear selection state on the sample set summary
      * drop index & clear selection state of each sample in the sampleset

    Delete a sample from a sample set, will:
      * drop index & clear selection state on the sample set summary
      * drop index & clear selection state on the sample

    Add a sample to a sample set, will:
      * drop index & clear selection state on the sample set summary
      * will NOT touch the selection state for the sample

Note that crossfilter indices are lazy created, as needed.
*/

export const genesetDelete = (genesetName) => (dispatch, getState) => {
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
    type: "sampleset: delete",
    genesetName,
    obsCrossfilter,
    annoMatrix: obsCrossfilter.annoMatrix,
  });
};

export const genesetAddGenes =
  (genesetName, genes) => async (dispatch, getState) => {
    const state = getState();
    const { obsCrossfilter: prevObsCrossfilter, annoMatrix } = state;
    const { schema } = annoMatrix;
    const varIndex = schema.annotations.var.index;
    const df = await annoMatrix.fetch("var", varIndex);
    const geneNames = df.col(varIndex).asArray();
    genes = genes.reduce((acc, gene) => {
      if (geneNames.indexOf(gene.geneSymbol) === -1) {
        postUserErrorToast(
          `${gene.geneSymbol} doesn't appear to be a valid sample name.`
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
      type: "sampleset: add samples",
      genesetName,
      genes,
      obsCrossfilter,
      annoMatrix: obsCrossfilter.annoMatrix,
    });
  };

export const genesetDeleteGenes =
  (genesetName, geneSymbols) => (dispatch, getState) => {
    const state = getState();
    const obsCrossfilter = dropGeneset(
      dispatch,
      state,
      genesetName,
      geneSymbols
    );
    return dispatch({
      type: "sampleset: delete samples",
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
