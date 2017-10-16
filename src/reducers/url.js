import uri from "urijs";

/************************************************************************
*************************************************************************
# This is all of the state that will be stored in the URL query params.
# It is reasonably important it be kept in the same place,
because if initial load or back button etc, we'll hear about that
and manually reconstruct the state from the url, overriding whatever
we had. We could of course listen for the "url changed" type in another
reducer, if this gets messy.
*************************************************************************
************************************************************************/

const url = (
  state = {},
  action
) => {
  switch (action.type) {
    case "url changed": {
      const {url} = action;

      // const parsed = uri.parseQuery(window.location.search);

      // console.log("onURLchange sees: ", state, parsed)

      return state;
    }
    case "category changed": {
      /* TODO catch categories here and store a representation of them */
      return state;
    }
    default:
      return state;
  }

  return state;
}

export default url;
