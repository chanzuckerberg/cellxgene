import * as globals from "../globals";
import URI from "urijs";

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
const attemptCategoricalMetadataSelection = (metadataField, value) => {
  return (dispatch, getState) => {
    dispatch({type: "categorical metadata filter selected start"})

    let uri = new URI()
    uri.setSearch(getState().selectedMetadata)
    uri.addSearch({[metadataField]: [value]})

    dispatch(
      requestCells(uri.search())
    ).then((res) => {
      if (res.error) {
        dispatch({type: "categorical metadata filter selected error"})
      } else {
        dispatch({type: "categorical metadata filter selected success", metadataField, value})
      }
    })
  }
}

/* DESELECT */
const attemptCategoricalMetadataDeselection = (metadataField, value) => {
  return (dispatch, getState) => {
    dispatch({type: "categorical metadata filter deselected start"})

    let uri = new URI()
    uri.setSearch(getState().selectedMetadata)
    uri.removeSearch({[metadataField]: [value]})

    dispatch(requestCells(uri.search())).then((res) => {
      if (res.error) {
        dispatch({type: "categorical metadata filter deselected error"})
      } else {
        dispatch({type: "categorical metadata filter deselected success", metadataField, value})
      }
    })
  }
}

const initialize = () => {
  return (dispatch, getState) => {
    dispatch({type: "initialize started"})
    fetch(`${globals.API.prefix}${globals.API.version}cells`, {
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

const ___hardcoded___requestGeneExpressionCountsPOST = () => {
  return (dispatch, getState) => {
    dispatch({type: "get expression started"})
    fetch(`${globals.API.prefix}${globals.API.version}expression`, {
        method: "POST",
        body: JSON.stringify({
          "genelist": [
            "AAK1",
            "1/2-SBSRNA4",
          ]
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

const requestDifferentialExpression = (celllist1, celllist2, num_genes = 5) => {
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
        data => dispatch({type: "request differential expression success", data}),
        error => dispatch({type: "request differential expression error", error})
      )
  }
}

export default {
  initialize,
  requestCells,
  requestGeneExpressionCounts,
  ___hardcoded___requestGeneExpressionCountsPOST,
  attemptCategoricalMetadataSelection,
  attemptCategoricalMetadataDeselection,
  requestDifferentialExpression
}
