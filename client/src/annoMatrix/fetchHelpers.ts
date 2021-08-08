export { doBinaryRequest, doFetch } from "../util/actionHelpers";

/* double URI encode - needed for query-param filters */
export function _dubEncURIComp(s: string | number | boolean): string {
  return encodeURIComponent(encodeURIComponent(s));
}

/* currently unused, consider deleting */
// eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any -- - FIXME: disabled temporarily on migrate to TS.
export function _fetchResult(promise: any) {
  let _status = "pending";
  const res = promise.then(
    // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
    (r: any) => {
      _status = "success";
      return r;
    },
    // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
    (e: any) => {
      _status = "error";
      throw e;
    }
  );

  res.status = () => _status;

  return res;
}
