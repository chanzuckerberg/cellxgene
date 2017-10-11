const updateURLMiddleware = store => next => action => {
  const oldState = store.getState();
  const nextAction = next(action);

  console.log('middleware', store, next, action)

  if (action.type === 'url changed') {
    return nextAction;
  }

  const state = store.getState();

  // Internal helper for working with URIs
  const oldURI = new URI(window.location.href);
  const newURI = new URI(oldURI).setQueryData({});

  newURI.setPath('/foo/bar');

  // Set the path based on state
  if (!state.isOnLandingPage && state.project.id) {
    newURI.setPath(newURI.getPath() + state.project.id + '/');
    newURI.addQueryData('baz', state.mode);
    newURI.addQueryData('bat', state.selection.activePageID);
  } else {
    newURI.setPath(newURI.getPath() + state.landingSection + '/');
  }

  // Avoid URL thrashing by replacing state while loading instead of pushing
  const newPath = newURI.toString();
  const oldPath = oldURI.toString();
  if (newPath !== oldPath) {
    if (
      (oldState.mode === 'asdf' &&
        state.mode === 'asdf' &&
        !oldState.isOnLandingPage) ||
      oldState.isLoadingProject !== state.isLoadingProject
    ) {
      window.history.replaceState(null, null, newPath);
    } else {
      window.history.pushState(null, null, newPath);
    }
  }

  return nextAction;
};

export default updateURLMiddleware;
