/*
action creators related to embeddings choice
*/

import { Action, ActionCreator } from "redux";
import { ThunkAction } from "redux-thunk";
import { AnnoMatrixObsCrossfilter } from "../annoMatrix";
import type { AppDispatch, RootState } from "../reducers";
import { _setEmbeddingSubset } from "../util/stateManager/viewStackHelpers";
import { Field } from "../common/types/schema";

// eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types --- FIXME: disabled temporarily on migrate to TS.
export async function _switchEmbedding(
  // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any -- - FIXME: disabled temporarily on migrate to TS.
  prevAnnoMatrix: any,
  // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any -- - FIXME: disabled temporarily on migrate to TS.
  prevCrossfilter: any,
  // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any -- - FIXME: disabled temporarily on migrate to TS.
  newEmbeddingName: any
) {
  /*
  DRY helper used by embedding action creators
  */
  const base = prevAnnoMatrix.base();
  const embeddingDf = await base.fetch("emb", newEmbeddingName);
  const annoMatrix = _setEmbeddingSubset(prevAnnoMatrix, embeddingDf);
  const obsCrossfilter = await new AnnoMatrixObsCrossfilter(
    annoMatrix,
    prevCrossfilter.obsCrossfilter
  ).select(Field.emb, newEmbeddingName, {
    mode: "all",
  });
  return [annoMatrix, obsCrossfilter];
}

export const layoutChoiceAction: ActionCreator<
  ThunkAction<Promise<void>, RootState, never, Action<"set layout choice">>
> =
  (newLayoutChoice: string) =>
  async (dispatch: AppDispatch, getState: () => RootState): Promise<void> => {
    /*
  On layout choice, make sure we have selected all on the previous layout, AND the new
  layout.
  */
    const { annoMatrix: prevAnnoMatrix, obsCrossfilter: prevCrossfilter } =
      getState();
    const [annoMatrix, obsCrossfilter] = await _switchEmbedding(
      prevAnnoMatrix,
      prevCrossfilter,
      newLayoutChoice
    );
    dispatch({
      type: "set layout choice",
      layoutChoice: newLayoutChoice,
      obsCrossfilter,
      annoMatrix,
    });
  };
