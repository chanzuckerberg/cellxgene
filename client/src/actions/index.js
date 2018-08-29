// jshint esversion: 6
import URI from "urijs";
import _ from "lodash";
import memoize from "memoize-one";
import * as globals from "../globals";
import store from "../reducers";

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

async function doRequestInitialize() {
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

async function doRequestCells(query) {
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

function doInitialDataLoad(query = "") {
  return catchErrorsWrap(async (dispatch, getState) => {
    dispatch({ type: "initial data load start" });
    dispatch({ type: "initialize started" });
    const rqstInit = doRequestInitialize();
    dispatch({ type: "request cells started" });
    const rqstCells = doRequestCells(query);

    let thereWasAnError = false;

    try {
      const data = await rqstInit;
      dispatch({ type: "initialize success", data });
    } catch (error) {
      thereWasAnError = true;
      dispatch({ type: "initialize error", error });
      dispatch({ type: "initial data load error", error });
    }

    if (!thereWasAnError) {
      try {
        const data = await rqstCells;
        dispatch({ type: "request cells success", data });
      } catch (error) {
        thereWasAnError = true;
        dispatch({ type: "request cells error", error });
        dispatch({ type: "initial data load error", error });
      }
    }

    if (!thereWasAnError) {
      dispatch({
        type: "initial data load complete (universe exists)",
        universe: getState().universe
      });
    }
  });
}

// /* SELECT */
// /* XXX TODO - this has not been refacotred for new redux state and is known to be broken */
// const regraph = () => {
//   return (dispatch, getState) => {
//     dispatch({ type: "regraph started" });
//
//     const state = getState();
//     const selectedMetadata = {};
//
//     _.each(state.controls2.categoricalAsBooleansMap, (options, field) => {
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

const regraph = () => dispatch =>
  dispatch({ type: "set World to current selection" });

const resetGraph = () => (dispatch, getState) =>
  dispatch({
    type: "reset World to eq Universe",
    universe: getState().universe
  });

// This code defends against the case where /expression returns a cellname
// never seen before (ie, not returned by /cells).   This should not happen
// (see https://github.com/chanzuckerberg/cellxgene-rest-api/issues/34) but
// occasionally does.
//
// XXX TODO - this code is only relevant in v0.1 REST API, and can be retired
// when we port to 0.2.
//
var makeMetadataMap = memoize(metadata => _.keyBy(metadata, "CellName"));
function cleanupExpressionResponse(data) {
  const s = store.getState();
  const metadata = makeMetadataMap(s.universe.obsAnnotations);
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
  const genesToFetch = _.filter(
    genes,
    g => !state.controls2.world.varDataByName(g)
  );

  if (!genesToFetch.length) {
    return dispatch({ type: "expression load not needed, all are cached" });
  }

  dispatch({ type: "expression load start" });
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
    return dispatch({ type: "expression load success", data });
  } catch (error) {
    dispatch({ type: "expression load error", error });
    throw error; // rethrow
  }
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
          [gene]: getState().controls2.world.varDataByName(gene)
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
          res[gene] = getState().controls2.world.varDataByName(gene);
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

export default {
  regraph,
  resetGraph,
  requestSingleGeneExpressionCountsForColoringPOST,
  requestDifferentialExpression,
  doInitialDataLoad
};
