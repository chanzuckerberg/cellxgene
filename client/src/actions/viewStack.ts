/*
The following actions manage the view stack for annoMatrix.

Conventions used and assumed elsewhere in the code base:
  * there will be zero or one clip view, and it will be the TOP view always.
  * there will be zero or more subset views

In other words, in our current use, we do not stack multiple clip views but we do
stack multiple subsets.

If these conventions change, code elsewhere (eg. menubar/clip.js) will need to
change as well.
*/
import { AnnoMatrixObsCrossfilter } from "../annoMatrix";
import {
  _clipAnnoMatrix,
  _userSubsetAnnoMatrix,
  _userResetSubsetAnnoMatrix,
} from "../util/stateManager/viewStackHelpers";

export const clipAction = (min: any, max: any) => (
  dispatch: any,
  getState: any
) => {
  /*
  apply a clip to the current annoMatrix.  By convention, the clip
  view is ALWAYS the top view.
  */
  const { annoMatrix: prevAnnoMatrix } = getState();
  const annoMatrix = _clipAnnoMatrix(prevAnnoMatrix, min, max);
  const obsCrossfilter = new AnnoMatrixObsCrossfilter(annoMatrix);
  dispatch({
    type: "set clip quantiles",
    clipQuantiles: { min, max },
    annoMatrix,
    obsCrossfilter,
  });
};

export const subsetAction = () => (dispatch: any, getState: any) => {
  /*
  Subset the annoMatrix to the current crossfilter selection by pushing a 
  subset view.

  By convention, a clip view is ALWAYS the top view, so if present, pop
  off and re-apply
  */
  const {
    annoMatrix: prevAnnoMatrix,
    obsCrossfilter: prevObsCrossfilter,
  } = getState();
  const annoMatrix = _userSubsetAnnoMatrix(
    prevAnnoMatrix,
    prevObsCrossfilter.allSelectedMask()
  );
  const obsCrossfilter = new AnnoMatrixObsCrossfilter(annoMatrix);
  dispatch({
    type: "subset to selection",
    annoMatrix,
    obsCrossfilter,
  });
};

export const resetSubsetAction = () => (dispatch: any, getState: any) => {
  /*
  Reset the annoMatrix to all data.  Because we may have multiple views
  stacked, we pop them all.  By convention, any clip transformation will
  be the top of the stack, and must be preserved.
  */

  const { annoMatrix: prevAnnoMatrix } = getState();
  const annoMatrix = _userResetSubsetAnnoMatrix(prevAnnoMatrix);
  const obsCrossfilter = new AnnoMatrixObsCrossfilter(annoMatrix);
  dispatch({
    type: "reset subset",
    annoMatrix,
    obsCrossfilter,
  });
};
