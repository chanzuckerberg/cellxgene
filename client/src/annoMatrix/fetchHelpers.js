export { doBinaryRequest } from "../util/actionHelpers";

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

  res.status = () => {
    return _status;
  };

  return res;
}

/* double URI encode - needed for query-param filters */
export function _dubEncURIComp(s) {
  return encodeURIComponent(encodeURIComponent(s));
}
