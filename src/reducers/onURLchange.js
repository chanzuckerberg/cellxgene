import URI from "urijs";
/*
  https://www.npmjs.com/package/url-parse
*/

const onURLchange = (
  state = {},
  action
) => {
  switch (action.type) {
    case 'url changed': {
      const {url} = action;

      const parsed = URI.parseQuery(window.location.search);

      console.log("onURLchange sees: ", parsed)

      return state;
    }
  }

  return state;
}

export default onURLchange;
