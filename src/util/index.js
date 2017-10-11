import URI from "urijs";

export const updateCategoricalQueryParams = (category, value) => {

  let newURL;

  /*
    there are four cases
    1. there are no query params
    2. there are query params, but our category isn't there yet
    3. there are query params and this is not the first entry to the given category (urijs may handle this with 'normalize')
    4. there are query params and ours already exists (remove it)
  */

  const uri = new URI(window.location.href)

  if (window.location.search === "") {
    newURL = uri.addQuery(category, value).toString(); /* #1 */
  } else if (uri.hasQuery(category, value) || uri.hasQuery(category, value, true)) { /* true param here means check arrays as well http://medialize.github.io/URI.js/docs.html#search-has */
    newURL = uri.removeQuery(category, value).toString(); /* #4 */
  } else {
    newURL = uri.addQuery(category, value).toString(); /* #2 & #3 are handled by URI */
  }

  window.history.pushState("", "", newURL)

}

export const updateContinuousParams = (category, range) => {
  /*
    this is similar to updateCategoricalQueryParams, but has a different check since
    in this case, same key (Unique_reads) different values (>500 vs >5000) should always replace previous range.
  */
}
