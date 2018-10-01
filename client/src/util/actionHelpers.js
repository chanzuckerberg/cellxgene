/*
Catch unexpected errors and make sure we don't lose them!
*/
export function catchErrorsWrap(fn) {
  return (dispatch, getState) => {
    fn(dispatch, getState).catch(error => {
      console.error(error);
      dispatch({ type: "UNEXPECTED ERROR", error });
    });
  };
}

/*
Bootstrap application with the initial data loading.
  * /config - application configuration
  * /schema - schema of dataframe
  * /annotations/obs - all metadata annotation
*/
export const doJsonRequest = async url => {
  const res = await fetch(url, {
    method: "get",
    headers: new Headers({
      "Content-Type": "application/json",
      "Accept-Encoding": "gzip, deflate, br"
    })
  });
  return res.json();
};
