// jshint esversion: 6
import _ from "lodash";
import memoize from "memoize-one";
import * as globals from "../globals";
import store from "../reducers";
import { Universe } from "../util/stateManager";

/*
Catch unexpected errors and make sure we don't lose them!
*/
function catchErrorsWrap(fn) {
  return (dispatch, getState) => {
    fn(dispatch, getState).catch(error => {
      console.error(error);
      dispatch({ type: "UNEXPECTED ERROR", error });
    });
  };
}

async function doRequestInitializeV01() {
  const res = await fetch(
    `${globals.API.prefix}${globals.API.version}initialize`,
    {
      method: "get",
      headers: new Headers({
        "Content-Type": "application/json"
      })
    }
  );
  return res.json();
}

async function doRequestCellsV01(query) {
  const res = await fetch(
    `${globals.API.prefix}${globals.API.version}cells${query}`,
    {
      method: "get",
      headers: new Headers({
        "Content-Type": "application/json"
      })
    }
  );
  return res.json();
}

function doInitialDataLoadV01(query = "") {
  return catchErrorsWrap(async dispatch => {
    dispatch({ type: "initial data load start" });
    try {
      const res = await Promise.all([
        doRequestInitializeV01(),
        doRequestCellsV01(query)
      ]);
      const universe = Universe.createUniverseFromRESTv01Response(
        res[0],
        res[1]
      );
      dispatch({
        type: "initial data load complete (universe exists)",
        universe
      });
    } catch (error) {
      dispatch({ type: "initial data load error", error });
    }
  });
}

/*
Bootstrap application with the initial data loading.
  * /config - application configuration
  * /schema - schema of dataframe
  * /annotations/obs - all metadata annotation
*/
const doJsonRequest = async url => {
  const res = await fetch(url, {
    method: "get",
    headers: new Headers({
      "Content-Type": "application/json"
    })
  });
  return res.json();
};

const doInitialDataLoadRESTV02 = () =>
  catchErrorsWrap(async dispatch => {
    dispatch({ type: "initial data load start" });

    try {
      const requests = _([
        "config",
        "schema",
        "annotations/obs",
        // TODO: "annotations/var",
        "layout/obs"
      ])
        .map(r => `${globals.API.prefix}${globals.API.version}${r}`)
        .map(url => doJsonRequest(url))
        .value();
      const results = await Promise.all(requests);
      const universe = Universe.createUniverseFromRestV02Response(...results);
      dispatch({
        type: "configuration load complete",
        config: requests[0]
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

// This code defends against the case where /expression returns a cellname
// never seen before (ie, not returned by /cells).   This should not happen
// (see https://github.com/chanzuckerberg/cellxgene-rest-api/issues/34) but
// occasionally does.
//
// XXX TODO - this code is only relevant in v0.1 REST API, and can be retired
// when we port to 0.2.
//
const makeMetadataMap = memoize(metadata => _.keyBy(metadata, "CellName"));
function cleanupExpressionResponse(data) {
  const s = store.getState();
  const { universe } = s.controls;
  const metadata = makeMetadataMap(universe.obsAnnotations);
  let errorFound = false;
  data.data.cells = _.filter(data.data.cells, cell => {
    if (!errorFound && !metadata[cell.cellname]) {
      errorFound = true;
      console.error(
        "Warning: /expression REST API returned unexpected cell names -- discarding surprises."
      );
    }
    return metadata[cell.cellname];
  });

  return data;
}

/*
Fetch [gene, ...] from V0.1 API.  Not an action function - just a helper
which implements the new expression data caching.
*/
async function _doRequestExpressionData(dispatch, getState, genes) {
  const state = getState();
  /* check cache and only fetch data we do not already have */
  const { universe } = state.controls;
  const genesToFetch = _.filter(genes, g => !universe.varDataCache[g]);

  dispatch({ type: "expression load start" });
  let expressionData = {}; // { gene: data }
  if (genesToFetch.length) {
    try {
      const res = await fetch(
        `${globals.API.prefix}${globals.API.version}expression`,
        {
          method: "POST",
          body: JSON.stringify({
            genelist: genes
          }),
          headers: new Headers({
            accept: "application/json",
            "Content-Type": "application/json"
          })
        }
      );
      let data = await res.json();
      data = cleanupExpressionResponse(data);
      data = Universe.convertExpressionRESTv01ToObject(universe, data);
      expressionData = {
        ...expressionData,
        ...data
      };
    } catch (error) {
      dispatch({ type: "expression load error", error });
      throw error; // rethrow
    }
  }

  // add the cached values
  _.forEach(genes, g => {
    if (expressionData[g] === undefined) {
      expressionData[g] = universe.varDataCache[g];
    }
  });

  return dispatch({ type: "expression load success", expressionData });
}

function requestSingleGeneExpressionCountsForColoringPOST(gene) {
  return async (dispatch, getState) => {
    dispatch({ type: "get single gene expression for coloring started" });
    try {
      await _doRequestExpressionData(dispatch, getState, [gene]);
      dispatch({
        type: "color by expression",
        gene,
        data: {
          [gene]: getState().controls.world.varDataCache[gene]
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
    await _doRequestExpressionData(dispatch, getState, genes);
    return dispatch({
      type: "get expression success",
      genes,
      data: _.transform(
        genes,
        (res, gene) => {
          res[gene] = getState().controls.world.varDataCache[gene];
        },
        {}
      )
    });
  } catch (error) {
    return dispatch({ type: "get expression error", error });
  }
};

const requestDifferentialExpression = (
  celllist1,
  celllist2,
  num_genes = 7
) => dispatch => {
  dispatch({ type: "request differential expression started" });
  fetch(`${globals.API.prefix}${globals.API.version}diffexpression`, {
    method: "POST",
    body: JSON.stringify({
      celllist1,
      celllist2,
      num_genes
    }),
    headers: new Headers({
      accept: "application/json",
      "Content-Type": "application/json"
    })
  })
    .then(res => res.json())
    .then(
      data => {
        /*
          kick off a secondary action to get all expression counts for all cells
          now that we know what the top expressed are
          */
        dispatch(
          requestGeneExpressionCountsPOST(
            _.union(data.data.celllist1.topgenes, data.data.celllist2.topgenes)
          )
        );
        /* then send the success case action through */
        return dispatch({
          type: "request differential expression success",
          data
        });
      },
      error =>
        dispatch({
          type: "request differential expression error",
          error
        })
    );
};

const doInitialDataLoad = doInitialDataLoadRESTV02;
export default {
  regraph,
  resetGraph,
  requestSingleGeneExpressionCountsForColoringPOST,
  requestDifferentialExpression,
  doInitialDataLoad
};
