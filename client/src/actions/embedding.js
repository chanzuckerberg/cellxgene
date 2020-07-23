/*
action creators related to embeddings choice
*/

import { AnnoMatrixObsCrossfilter } from "../annoMatrix";
import { _setEmbeddingSubset } from "../util/stateManager/viewStackHelpers";

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
    obsCrossfilter: prevObsCrossfilter,
    layoutChoice,
  } = getState();

  const embeddingDf = await prevAnnoMatrix.base().fetch("emb", newLayoutChoice);
  const annoMatrix = _setEmbeddingSubset(prevAnnoMatrix, embeddingDf);
  const obsCrossfilter = await new AnnoMatrixObsCrossfilter(annoMatrix).select(
    "emb",
    newLayoutChoice,
    {
      mode: "all",
    }
  );
  dispatch({
    type: "set layout choice",
    layoutChoice: newLayoutChoice,
    obsCrossfilter,
    annoMatrix,
  });
};
