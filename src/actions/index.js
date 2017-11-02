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

export default {
  initialize,
  requestCells,
  attemptCategoricalMetadataSelection,
  attemptCategoricalMetadataDeselection
}
