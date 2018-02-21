import * as globals from "../globals";
import URI from "urijs";
import _ from "lodash";

const requestCells = (query = "") => {
  return (dispatch) => {
    dispatch({type: "request cells started"})
    return fetch(`${globals.API.prefix}${globals.API.version}cells${query}`, {
        method: "get",
        headers: new Headers({
          'Content-Type': 'application/json'
        })
      })
      .then(res => res.json())
      .then(
        data => dispatch({type: "request cells success", data}),
        error => dispatch({type: "request cells error", error})
      )
  }
}

/* SELECT */
const regraph = () => {
  return (dispatch, getState) => {
    dispatch({type: "regraph started"})

    const state = getState()
    const selectedMetadata = {};

    _.each(state.controls.categoricalAsBooleansMap, (options, field) => {
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
              selectedMetadata[field].push(option)
            } else if (!selectedMetadata[field]) {
              selectedMetadata[field] = [option]
            }
          }
        })
      }
    })

    let uri = new URI()
    uri.setSearch(selectedMetadata)
    console.log(uri.search(), selectedMetadata)

    dispatch(
      requestCells(uri.search())
    ).then((res) => {
      if (res.error) {
        dispatch({type: "regraph success"})
      } else {
        dispatch({type: "regraph error"})
      }
    })
  }
}

const initialize = () => {
  return (dispatch, getState) => {
    dispatch({type: "initialize started"})
    fetch(`${globals.API.prefix}${globals.API.version}initialize`, {
        method: "get",
        headers: new Headers({
          'Content-Type': 'application/json'
        })
      })
      .then(res => res.json())
      .then(
        data => dispatch({type: "initialize success", data}),
        error => dispatch({type: "initialize error", error})
      )

  }
}

const requestGeneExpressionCounts = () => {
  return (dispatch, getState) => {
    dispatch({type: "get expression started"})
    fetch(`${globals.API.prefix}${globals.API.version}expression`, {
        method: "get",
        headers: new Headers({
          'accept': 'application/json'
        })
      })
      .then(res => res.json())
      .then(
        data => dispatch({type: "get expression success", data}),
        error => dispatch({type: "get expression error", error})
      )
  }
}

const requestSingleGeneExpressionCountsForColoringPOST = (gene) => {
  return (dispatch, getState) => {
    dispatch({type: "get single gene expression for coloring started"})
    fetch(`${globals.API.prefix}${globals.API.version}expression`, {
        method: "POST",
        body: JSON.stringify({
          "genelist": [gene]
        }),
        headers: new Headers({
          "accept": "application/json",
          "Content-Type": "application/json"
        })
      })
      .then(res => res.json())
      .then(
        data => dispatch({
          type: "color by expression",
          gene: gene,
          data
        }),
        error => dispatch({type: "get single gene expression for coloring error", error})
      )
  }
}

const requestGeneExpressionCountsPOST = (genes) => {
  return (dispatch, getState) => {
    dispatch({type: "get expression started"})
    fetch(`${globals.API.prefix}${globals.API.version}expression`, {
        method: "POST",
        body: JSON.stringify({
          "genelist": genes
        }),
        headers: new Headers({
          "accept": "application/json",
          "Content-Type": "application/json"
        })
      })
      .then(res => res.json())
      .then(
        data => dispatch({type: "get expression success", data}),
        error => dispatch({type: "get expression error", error})
      )
  }
}

const requestDifferentialExpression = (celllist1, celllist2, num_genes = 7) => {
  return (dispatch, getState) => {
    dispatch({type: "request differential expression started"})
    fetch(`${globals.API.prefix}${globals.API.version}diffexpression`, {
        method: "POST",
        body: JSON.stringify({
          celllist1,
          celllist2,
          num_genes,
        }),
        headers: new Headers({
          "accept": "application/json",
          "Content-Type": "application/json"
        })
      })
      .then(res => res.json())
      .then(
        (data) => {
          /* kick off a secondary action to get all expression counts for all cells now that we know what the top expressed are */
          dispatch(
            requestGeneExpressionCountsPOST(
              _.union(data.data.celllist1.topgenes, data.data.celllist2.topgenes) // ["GPM6B", "FEZ1", "TSPAN31", "PCSK1N", "TUBA1A", "GPM6A", "CLU", "FCER1G", "TYROBP", "C1QB", "CD74", "CYBA", "GPX1", "TMSB4X"]
            )
          )
          /* then send the success case action through */
          return dispatch({type: "request differential expression success", data})
        },
        error => dispatch({type: "request differential expression error", error})
      )
  }
}

export default {
  initialize,
  requestCells,
  regraph,
  requestGeneExpressionCounts,
  requestGeneExpressionCountsPOST,
  requestSingleGeneExpressionCountsForColoringPOST,
  requestDifferentialExpression
}
