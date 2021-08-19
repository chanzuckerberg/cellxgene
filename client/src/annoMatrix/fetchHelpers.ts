export { doBinaryRequest, doFetch } from "../util/actionHelpers";

/* double URI encode - needed for query-param filters */
export function _dubEncURIComp(s: any) {
  return encodeURIComponent(encodeURIComponent(s));
}

/* currently unused, consider deleting */
export function _fetchResult(promise: any) {
  let _status = "pending";
  const res = promise.then(
    (r: any) => {
      _status = "success";
      return r;
    },
    (e: any) => {
      _status = "error";
      throw e;
    }
  );

  res.status = () => {
    return _status;
  };

  return res;
}
