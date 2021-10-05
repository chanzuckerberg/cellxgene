/*
action creators related to embeddings choice
*/

import { AnnoMatrixObsCrossfilter } from "../annoMatrix";
import { _setEmbeddingSubset } from "../util/stateManager/viewStackHelpers";
import { API } from "../globals";

export async function _switchEmbedding(
  prevAnnoMatrix,
  prevCrossfilter,
  oldEmbeddingName,
  newEmbeddingName
) {
  /*
  DRY helper used by this and reembedding action creators
  */
  const base = prevAnnoMatrix.base();
  const embeddingDf = await base.fetch("emb", newEmbeddingName);

  const annoMatrix = _setEmbeddingSubset(prevAnnoMatrix, embeddingDf);
  let obsCrossfilter = await new AnnoMatrixObsCrossfilter(
    annoMatrix,
    prevCrossfilter.obsCrossfilter
  ).dropDimension("emb", oldEmbeddingName);
  obsCrossfilter = await obsCrossfilter.dropDimension("obs", "name_0");
  obsCrossfilter = await obsCrossfilter.select("emb", newEmbeddingName, {
    mode: "all",
  });
  return [annoMatrix, obsCrossfilter];
}

export const requestDeleteEmbedding = (toDelete) => async (
    _dispatch,
    _getState
  ) => {
  const res = await fetch(
    `${API.prefix}${API.version}layout/obsm`,
    {
      method: "PUT",
      headers: new Headers({
        Accept: "application/octet-stream",
        "Content-Type": "application/json",
      }),
      body: JSON.stringify({
        embNames: toDelete,
      }),
      credentials: "include",
    }
  );
  const schema = await res.json();

  if (res.ok && res.headers.get("Content-Type").includes("application/json")) {      
    return schema;
  }

  // else an error
  let msg = `Unexpected HTTP response ${res.status}, ${res.statusText}`;
  const body = await res.text();
  if (body && body.length > 0) {
    msg = `${msg} -- ${body}`;
  }
  throw new Error(msg);
}

export const layoutChoiceAction = (newLayoutChoice) => async (
  dispatch,
  getState
) => {
  /*
  On layout choice, make sure we have selected all on the previous layout, AND the new
  layout.
  */
  const {
    annoMatrix: prevAnnoMatrix,
    obsCrossfilter: prevCrossfilter,
    layoutChoice,
  } = getState();
  const [annoMatrix, obsCrossfilter] = await _switchEmbedding(
    prevAnnoMatrix,
    prevCrossfilter,
    layoutChoice.current,
    newLayoutChoice
  );
  dispatch({
    type: "set layout choice",
    layoutChoice: newLayoutChoice,
    obsCrossfilter,
    annoMatrix,
  });
};
