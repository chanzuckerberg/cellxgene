// jshint esversion: 6
import _ from "lodash";
import * as globals from "../globals";
import { Universe, kvCache } from "../util/stateManager";
import {
  catchErrorsWrap,
  doJsonRequest,
  rangeEncodeIndices
} from "../util/actionHelpers";

/*
Bootstrap application with the initial data loading.
  * /config - application configuration
  * /schema - schema of dataframe
  * /annotations/obs - all metadata annotation
*/
const doInitialDataLoad = () =>
  catchErrorsWrap(async dispatch => {
    dispatch({ type: "initial data load start" });

    try {
      const requests = _([
        "config",
        "schema",
        "annotations/obs",
        "annotations/var",
        "layout/obs"
      ])
        .map(r => `${globals.API.prefix}${globals.API.version}${r}`)
        .map(url => doJsonRequest(url))
        .value();
      const results = await Promise.all(requests);

      /* set config defaults */
      const config = { ...globals.configDefaults, ...results[0].config };
      const [, schema, obsAnno, varAnno, obsLayout] = [...results];
      const universe = Universe.createUniverseFromRestV02Response(
        config,
        schema,
        obsAnno,
        varAnno,
        obsLayout
      );

      dispatch({
        type: "configuration load complete",
        config
      });
      dispatch({
        type: "initial data load complete (universe exists)",
        universe
      });
    } catch (error) {
      dispatch({ type: "initial data load error", error });
    }
  });

// XXX TODO - this is the old code for doing a regraph.  Preserving it solely
// until we port to 0.2 API.   The new UX for regraph can't be implemented on
// the 0.1 API (doesn't allow for re-layout on arbitrary sets of cells), so just
// punting for now.  See ticket #88
//
//
// /* SELECT */
// const regraph = () => {
//   return (dispatch, getState) => {
//     dispatch({ type: "regraph started" });
//
//     const state = getState();
//     const selectedMetadata = {};
//
//     _.each(state.controls.categoricalAsBooleansMap, (options, field) => {
//       let atLeastOneOptionDeselected = false;
//
//       _.each(options, (isActive, option) => {
//         if (!isActive) {
//           atLeastOneOptionDeselected = true;
//         }
//       });
//
//       if (atLeastOneOptionDeselected) {
//         _.each(options, (isActive, option) => {
//           if (isActive) {
//             if (selectedMetadata[field]) {
//               selectedMetadata[field].push(option);
//             } else if (!selectedMetadata[field]) {
//               selectedMetadata[field] = [option];
//             }
//           }
//         });
//       }
//     });
//
//     let uri = new URI();
//     uri.setSearch(selectedMetadata);
//     console.log(uri.search(), selectedMetadata);
//
//     dispatch(requestCells(uri.search())).then(res => {
//       if (res.error) {
//         dispatch({ type: "regraph error" });
//       } else {
//         dispatch({ type: "regraph success" });
//       }
//     });
//   };
// };

const regraph = () => (dispatch, getState) => {
  const { universe, world, crossfilter } = getState().controls;
  dispatch({
    type: "set World to current selection",
    universe,
    world,
    crossfilter
  });
};

const resetGraph = () => (dispatch, getState) =>
  dispatch({
    type: "reset World to eq Universe",
    universe: getState().controls.universe
  });

/*
Fetch expression vectors for each gene in genes.   This is NOT an action
function, but rather a helper to be called from an action helper that
needs expression data.

Transparently utilizes cached data if it is already present.
*/
async function _doRequestExpressionData(dispatch, getState, genes) {
  const state = getState();
  const { universe } = state.controls;
  /* preload data already in cache */
  let expressionData = _.transform(genes, (expData, g) => {
    const data = kvCache.get(universe.varDataCache, g);
    if (data) {
      expData[g] = data;
    }
  }); // --> { gene: data }
  /* make a list of genes for which we do not have data */
  const genesToFetch = _.filter(genes, g => expressionData[g] === undefined);

  dispatch({ type: "expression load start" });

  /* Fetch data for any genes not in cache */
  if (genesToFetch.length) {
    try {
      // XXX: TODO - this could be using /data/var rather than /data/obs,
      // as that would simplify the transformation in
      // convertExpressionRESTv02ToObject
      const res = await fetch(
        `${globals.API.prefix}${globals.API.version}data/obs`,
        {
          method: "PUT",
          body: JSON.stringify({
            filter: {
              var: {
                annotation_value: [{ name: "name", values: genesToFetch }]
              }
            }
          }),
          headers: new Headers({
            accept: "application/json",
            "Accept-Encoding": "gzip, deflate, br",
            "Content-Type": "application/json"
          })
        }
      );
      const data = await res.json();
      expressionData = {
        ...expressionData,
        ...Universe.convertExpressionRESTv02ToObject(universe, data)
      };
    } catch (error) {
      dispatch({ type: "expression load error", error });
      throw error; // rethrow
    }
  }

  dispatch({ type: "expression load success", expressionData });
  return expressionData;
}

function requestSingleGeneExpressionCountsForColoringPOST(gene) {
  return async (dispatch, getState) => {
    dispatch({ type: "get single gene expression for coloring started" });
    try {
      const expressionData = await _doRequestExpressionData(
        dispatch,
        getState,
        [gene]
      );
      dispatch({
        type: "color by expression",
        gene,
        data: {
          [gene]: expressionData[gene]
        }
      });
    } catch (error) {
      dispatch({
        type: "get single gene expression for coloring error",
        error
      });
    }
  };
}

const requestGeneExpressionCountsPOST = genes => async (dispatch, getState) => {
  dispatch({ type: "get expression started" });
  try {
    const expressionData = await _doRequestExpressionData(
      dispatch,
      getState,
      genes
    );
    return dispatch({
      type: "get expression success",
      genes,
      data: expressionData
    });
  } catch (error) {
    return dispatch({ type: "get expression error", error });
  }
};

const requestUserDefinedGene = gene => async dispatch => {
  dispatch({ type: "request user defined gene started" });
  try {
    const data = await dispatch(requestGeneExpressionCountsPOST([gene]));

    /* then send the success case action through */
    return dispatch({
      type: "request user defined gene success",
      data
    });
  } catch (error) {
    return dispatch({
      type: "request user defined gene error",
      error
    });
  }
};

const requestDifferentialExpression = (set1, set2, num_genes = 10) => async (
  dispatch,
  getState
) => {
  dispatch({ type: "request differential expression started" });
  try {
    /*
    Steps:
    1. get the most differentially expressed genes
    2. get expression data for each
    */
    const state = getState();
    const { universe } = state.controls;
    const set1ByIndex = rangeEncodeIndices(
      _.map(set1, s => universe.obsNameToIndexMap[s])
    );
    const set2ByIndex = rangeEncodeIndices(
      _.map(set2, s => universe.obsNameToIndexMap[s])
    );
    const diffExpFetch = await fetch(
      `${globals.API.prefix}${globals.API.version}diffexp/obs`,
      {
        method: "POST",
        headers: new Headers({
          Accept: "application/json",
          "Accept-Encoding": "gzip, deflate, br",
          "Content-Type": "application/json"
        }),
        body: JSON.stringify({
          mode: "topN",
          count: num_genes,
          set1: { filter: { obs: { index: set1ByIndex } } },
          set2: { filter: { obs: { index: set2ByIndex } } }
        })
      }
    );
    const data = await diffExpFetch.json();
    // result is [ [varIdx, ...], ... ]
    const topNGenes = _.map(data, r => universe.varAnnotations[r[0]].name);

    /*
    Kick off secondary action to fetch all of the expression data for the
    topN expressed genes.
    */
    await dispatch(requestGeneExpressionCountsPOST(topNGenes));

    /* then send the success case action through */
    return dispatch({
      type: "request differential expression success",
      data
    });
  } catch (error) {
    return dispatch({
      type: "request differential expression error",
      error
    });
  }
};

export default {
  regraph,
  resetGraph,
  requestSingleGeneExpressionCountsForColoringPOST,
  requestGeneExpressionCountsPOST,
  requestDifferentialExpression,
  requestUserDefinedGene,
  doInitialDataLoad
};
