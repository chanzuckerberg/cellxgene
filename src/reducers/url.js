import uri from "urijs";

/************************************************************************
*************************************************************************
1. This is all of the state that will be stored in the URL query params.
1a. It is reasonably important it be kept in the same place,
because if initial load or back button etc, we'll hear about that
and manually reconstruct the state from the url, overriding whatever
we had. We could of course listen for the "url changed" type in another
reducer, if this gets messy.
1b. The data structures used here are constrained as well, as urijs
needs to be able to parse them
*************************************************************************
************************************************************************/

const url = (
  state = {
    selectedMetadata: {}
  },
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

      return state;
    }
    default:
      return state;
  }

  return state;
}

export default url;
