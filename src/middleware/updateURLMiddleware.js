import uri from "urijs";

const updateURLMiddleware = (store) => {
  return (next) => {
    return (action) => {
      const oldState = store.getState();
      const nextAction = next(action);

      if (action.type === 'url changed') {
        /* we don't handle pop state here - we handle it in the url reducer */
        return nextAction;
      }

      const state = store.getState();



      // const oldURI = URI(window.location.href)
      // const newURI = URI(oldURI).setQuery({})

      // if (window.location.search === "") {
      //   newURL = uri.addQuery(category, value).toString(); /* #1 */
      // } else if (uri.hasQuery(category, value) || uri.hasQuery(category, value, true)) { /* true param here means check arrays as well http://medialize.github.io/URI.js/docs.html#search-has */
      //   newURL = uri.removeQuery(category, value).toString(); /* #4 */
      // } else {
      //   newURL = uri.addQuery(category, value).toString(); /* #2 & #3 are handled by URI */
      // }
      //
      // window.history.pushState("", "", newURL)


      //
      // // Internal helper for working with URIs
      // const oldURI = new URI(window.location.href);
      // const newURI = new URI(oldURI).setQueryData({});
      //
      // newURI.setPath('/foo/bar');
      //
      // // Set the path based on state
      // if (!state.isOnLandingPage && state.project.id) {
      //   newURI.setPath(newURI.getPath() + state.project.id + '/');
      //   newURI.addQueryData('baz', state.mode);
      //   newURI.addQueryData('bat', state.selection.activePageID);
      // } else {
      //   newURI.setPath(newURI.getPath() + state.landingSection + '/');
      // }
      //
      // // Avoid URL thrashing by replacing state while loading instead of pushing
      // const newPath = newURI.toString();
      // const oldPath = oldURI.toString();
      // if (newPath !== oldPath) {
      //   if (
      //     (oldState.mode === 'asdf' &&
      //       state.mode === 'asdf' &&
      //       !oldState.isOnLandingPage) ||
      //     oldState.isLoadingProject !== state.isLoadingProject
      //   ) {
      //     window.history.replaceState(null, null, newPath);
      //   } else {
      //     window.history.pushState(null, null, newPath);
      //   }
      // }

      return nextAction;
    }
  }
};

export default updateURLMiddleware;
