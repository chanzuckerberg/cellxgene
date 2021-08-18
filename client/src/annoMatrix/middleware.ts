/*
Garbage collection / cache management support

Middleware that knows how to pull annoMatrix from the undoable state,
and pass it along to the AnnoMatrix class for possible cache GC.

Private interface.

Future work item: this middleware knows internal details of both the
Undoable metareducer and the AnnoMatrix private API.  It would be helpful
to make the Undoable interface better factored.
*/

import { Action, Dispatch, MiddlewareAPI } from "redux";
import AnnoMatrix from "./annoMatrix";
import { GCHints } from "../common/types/entities";

const annoMatrixGC =
  (store: MiddlewareAPI) =>
  // GC middleware doesn't add any extra types to dispatch; it just executes GC and continues.
  (next: Dispatch) =>
  (action: Action): Action => {
    if (_itIsTimeForGC()) {
      _doGC(store);
    }
    return next(action);
  };

let lastGCTime = 0;
const InterGCDelayMS = 30 * 1000; // 30 seconds
function _itIsTimeForGC(): boolean {
  /*
  we don't want to run GC on every dispatch, so throttle it a bit.

  Runs every InterGCDelay period
  */
  const now = Date.now();
  if (now - lastGCTime > InterGCDelayMS) {
    lastGCTime = now;
    return true;
  }
  return false;
}

function _doGC(store: MiddlewareAPI): void {
  const state = store.getState();

  // these should probably be a function imported from undoable.js, etc, as
  // they have overly intimate knowledge of our reducers.
  const undoablePast = state["@@undoable/past"];
  const undoableFuture = state["@@undoable/future"];
  const undoableStack = undoablePast
    .concat(undoableFuture)
    // eslint-disable-next-line @typescript-eslint/no-explicit-any --- TODO revisit: waiting for typings from /reducers
    .flatMap((snapshot: any) =>
      snapshot
        // eslint-disable-next-line @typescript-eslint/no-explicit-any --- TODO revisit: waiting for typings from /reducers
        .filter((v: any) => v[0] === "annoMatrix")
        // eslint-disable-next-line @typescript-eslint/no-explicit-any --- TODO revisit: waiting for typings from /reducers
        .map((v: any) => v[1])
    );
  const currentAnnoMatrix = state.annoMatrix;

  /*
  We want to identify those matrixes currently "hot", ie, linked from the current annoMatrix,
  as our current gc algo is more aggressive with those not hot.
  */
  const allAnnoMatrices = new Map<AnnoMatrix, GCHints>(
    undoableStack.map((m: AnnoMatrix) => [m, { isHot: false }])
  );
  let am = currentAnnoMatrix;
  while (am && am.isView) {
    allAnnoMatrices.set(am, { isHot: true });
    am = am.viewOf;
  }
  allAnnoMatrices.forEach((hints, annoMatrix: AnnoMatrix) =>
    annoMatrix._gc(hints)
  );
}

export default annoMatrixGC;
