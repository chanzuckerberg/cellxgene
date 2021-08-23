export { doBinaryRequest, doFetch } from "../util/actionHelpers";

/* double URI encode - needed for query-param filters */
export function _dubEncURIComp(s) {
  return encodeURIComponent(encodeURIComponent(s));
}

/* currently unused, consider deleting */
export function _fetchResult(promise) {
  let _status = "pending";
  const res = promise.then(
    (r) => {
      _status = "success";
      return r;
    },
    (e) => {
      _status = "error";
      throw e;
    }
  );

  res.status = () => _status;

  return res;
}
