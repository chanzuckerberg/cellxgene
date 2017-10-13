import uri from "urijs";
/*
  https://www.npmjs.com/package/url-parse
*/

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
