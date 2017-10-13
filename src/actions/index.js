import * as globals from "../globals";

export const requestCells = (categoriesToSubsetBy) => {
  return (dispatch, getState) => {
    dispatch({type: "request cells started"})
    fetch(`${globals.API.prefix}${globals.API.version}cells`, {
        method: "get",
        headers: new Headers({
          'Content-Type': 'application/json'
        }),
        // body: JSON.stringify({
        //   // "celllist": ["1001000173.G8", "1001000173.D4"],
        //   "genelist": ["1/2-SBSRNA4", "A1BG", "A1BG-AS1", "A1CF", "A2LD1", "A2M", "A2ML1", "A2MP1", "A4GALT"]
        // })
      })
      .then(res => res.json())
      .then(
        data => dispatch({type: "request cells success", data}),
        error => dispatch({type: "request cells error", error})
      )

  }
}
