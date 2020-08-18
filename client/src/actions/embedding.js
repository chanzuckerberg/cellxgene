/*
action creators related to embeddings choice
*/

import { AnnoMatrixObsCrossfilter } from "../annoMatrix";
import { _setEmbeddingSubset } from "../util/stateManager/viewStackHelpers";

export async function _switchEmbedding(prevAnnoMatrix, newEmbeddingName) {
  /*
  DRY helper used by this and reembedding action creators
  */
  const base = prevAnnoMatrix.base();
  const embeddingDf = await base.fetch("emb", newEmbeddingName);
  const annoMatrix = _setEmbeddingSubset(prevAnnoMatrix, embeddingDf);
  const obsCrossfilter = await new AnnoMatrixObsCrossfilter(annoMatrix).select(
    "emb",
    newEmbeddingName,
    {
      mode: "all",
    }
  );
  return [annoMatrix, obsCrossfilter];
}

export const layoutChoiceAction = (newLayoutChoice) => async (
  dispatch,
  getState
) => {
  /*
  On layout choice, make sure we have selected all on the previous layout, AND the new
  layout.
  */
  const { annoMatrix: prevAnnoMatrix } = getState();
  const [annoMatrix, obsCrossfilter] = await _switchEmbedding(
    prevAnnoMatrix,
    newLayoutChoice
  );
  dispatch({
    type: "set layout choice",
    layoutChoice: newLayoutChoice,
    obsCrossfilter,
    annoMatrix,
  });
};
