// jshint esversion: 6
import * as globals from "../globals";
import store from "../reducers";
import URI from "urijs";
import _ from "lodash";
import memoize from "memoize-one";

/*
Catch unexpected errors and make sure we don't lose them!
*/
function catchErrors(fn) {
  return function(dispatch, getState) {
    fn(dispatch, getState).catch(error => {
      console.error(error);
      dispatch({ type: "UNEXPECTED ERROR", error });
    });
  };
}

async function doRequestInitialize() {
  let res = await fetch(
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
  let res = await fetch(
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
  return catchErrors(async function(dispatch, getState) {
    dispatch({ type: "initial data load start" });
    dispatch({ type: "initialize started" });
    let rqstInit = doRequestInitialize();
    dispatch({ type: "request cells started" });
    let rqstCells = doRequestCells(query);

    let thereWasAnError = false;

    try {
      let data = await rqstInit;
      dispatch({ type: "initialize success", data });
    } catch (error) {
      thereWasAnError = true;
      dispatch({ type: "initialize error", error });
      dispatch({ type: "initial data load error", error });
    }

    if (!thereWasAnError) {
      try {
        let data = await rqstCells;
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
        universe: getState().dataframe
      });
    }
  });
}

// XXX dead code
// const requestCells = (query = "") => {
//   return dispatch => {
//     dispatch({ type: "request cells started" });
//     return fetch(`${globals.API.prefix}${globals.API.version}cells${query}`, {
//       method: "get",
//       headers: new Headers({
//         "Content-Type": "application/json"
//       })
//     })
//       .then(res => res.json())
//       .then(
//         data => dispatch({ type: "request cells success", data }),
//         error => dispatch({ type: "request cells error", error })
//       );
//   };
// };

/* SELECT */
/* XXX TODO - this has not been refacotred for new redux state and is known to be broken */
const regraph = () => {
  return (dispatch, getState) => {
    dispatch({ type: "regraph started" });

    const state = getState();
    const selectedMetadata = {};

    _.each(state.controls2.categoricalAsBooleansMap, (options, field) => {
      let atLeastOneOptionDeselected = false;

      _.each(options, (isActive, option) => {
        if (!isActive) {
          atLeastOneOptionDeselected = true;
        }
      });

      if (atLeastOneOptionDeselected) {
        _.each(options, (isActive, option) => {
          if (isActive) {
            if (selectedMetadata[field]) {
              selectedMetadata[field].push(option);
            } else if (!selectedMetadata[field]) {
              selectedMetadata[field] = [option];
            }
          }
        });
      }
    });

    let uri = new URI();
    uri.setSearch(selectedMetadata);
    console.log(uri.search(), selectedMetadata);

    dispatch(requestCells(uri.search())).then(res => {
      if (res.error) {
        dispatch({ type: "regraph error" });
      } else {
        dispatch({ type: "regraph success" });
      }
    });
  };
};

const resetGraph = () => {
  return (dispatch, getState) => {
    dispatch({ type: "reset graph" });
  };
};

// XXX: dead code
// const initialize = () => {
//   return (dispatch, getState) => {
//     dispatch({ type: "initialize started" });
//     fetch(`${globals.API.prefix}${globals.API.version}initialize`, {
//       method: "get",
//       headers: new Headers({
//         "Content-Type": "application/json"
//       })
//     })
//       .then(res => res.json())
//       .then(
//         data => dispatch({ type: "initialize success", data }),
//         error => dispatch({ type: "initialize error", error })
//       );
//   };
// };

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
  const metadata = makeMetadataMap(s.dataframe.obsAnnotations);
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

// XXX dead code
// const requestGeneExpressionCounts = () => {
//   return (dispatch, getState) => {
//     dispatch({ type: "get expression started" });
//     fetch(`${globals.API.prefix}${globals.API.version}expression`, {
//       method: "get",
//       headers: new Headers({
//         accept: "application/json"
//       })
//     })
//       .then(res => res.json())
//       .then(data => cleanupExpressionResponse(data))
//       .then(
//         data => dispatch({ type: "get expression success", data }),
//         error => dispatch({ type: "get expression error", error })
//       );
//   };
// };

const requestSingleGeneExpressionCountsForColoringPOST = gene => {
  return (dispatch, getState) => {
    dispatch({ type: "get single gene expression for coloring started" });
    fetch(`${globals.API.prefix}${globals.API.version}expression`, {
      method: "POST",
      body: JSON.stringify({
        genelist: [gene]
      }),
      headers: new Headers({
        accept: "application/json",
        "Content-Type": "application/json"
      })
    })
      .then(res => res.json())
      .then(data => cleanupExpressionResponse(data))
      .then(
        data =>
          dispatch({
            type: "color by expression",
            gene: gene,
            data
          }),
        error =>
          dispatch({
            type: "get single gene expression for coloring error",
            error
          })
      );
  };
};

const requestGeneExpressionCountsPOST = genes => {
  return (dispatch, getState) => {
    dispatch({ type: "get expression started" });
    fetch(`${globals.API.prefix}${globals.API.version}expression`, {
      method: "POST",
      body: JSON.stringify({
        genelist: genes
      }),
      headers: new Headers({
        accept: "application/json",
        "Content-Type": "application/json"
      })
    })
      .then(res => res.json())
      .then(data => cleanupExpressionResponse(data))
      .then(
        data => dispatch({ type: "get expression success", data }),
        error => dispatch({ type: "get expression error", error })
      );
  };
};

const requestDifferentialExpression = (celllist1, celllist2, num_genes = 7) => {
  return (dispatch, getState) => {
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
          /* kick off a secondary action to get all expression counts for all cells now that we know what the top expressed are */
          dispatch(
            requestGeneExpressionCountsPOST(
              _.union(
                data.data.celllist1.topgenes,
                data.data.celllist2.topgenes
              ) // ["GPM6B", "FEZ1", "TSPAN31", "PCSK1N", "TUBA1A", "GPM6A", "CLU", "FCER1G", "TYROBP", "C1QB", "CD74", "CYBA", "GPX1", "TMSB4X"]
            )
          );
          /* then send the success case action through */
          return dispatch({
            type: "request differential expression success",
            data
          });
        },
        error =>
          dispatch({ type: "request differential expression error", error })
      );
  };
};

export default {
  // initialize,
  // requestCells,
  regraph,
  resetGraph,
  /* XXX: these are unused outside of this file */
  // requestGeneExpressionCounts,
  // requestGeneExpressionCountsPOST,
  requestSingleGeneExpressionCountsForColoringPOST,
  requestDifferentialExpression,
  doInitialDataLoad
};
