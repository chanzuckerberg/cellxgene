import queryString from "query-string";
import URI from "urijs";

export const updateCategoricalQueryParams = (category, value) => {

  let newURL;

  /*
    there are four cases
    1. there are no query params
    2. there are query params, but our category isn't there yet
    3. there are query params and this is not the first entry to the given category
    4. there are query params and ours already exists (remove it)
  */

  if (window.location.search === "") {
    /* #1 */
    newURL = `${window.location.href}?${queryString.stringify({category: [value]})}`
    window.history.pushState("", "", newUrl)
  } else {

    const existingParams = queryString.parse(window.location.search);
    const categoryAlreadyExists = existingParams[category];
    const valueAlreadyExists = existingParams[category] && categoryAlreadyExists[value];

    if (valueAlreadyExists) {
      /* #4 (removal) */

      newURL = `${window.location.href}?${queryString.stringify({category: [value]})}`
    } else {
      /* #3 */
      existingParams
    }

  }

}

export const updateContinuousParams = (category, range) => {

}
