export const createSession = (prompt, polygon) => (dispatch, getState) => {
  const { obsCrossfilter } = getState();
  const indexes = obsCrossfilter.allSelectedLabels();
  dispatch({
    type: "chat: create session",
    prompt,
    polygon,
    indexes,
  });
};
