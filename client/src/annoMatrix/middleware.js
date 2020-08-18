/*
Garbage collection / cache management support

Middleware that knows how to pull annoMatrix from the undoable state,
and pass it along to the AnnoMatrix class for possible cache GC.

Private interface.

Future work item: this middleware knows internal details of both the
Undoable metareducer and the AnnoMatrix private API.  It would be helpful
to make the Undoable interface better factored.
*/

const annoMatrixGC = (store) => (next) => (action) => {
  if (_itIsTimeForGC()) {
    _doGC(store);
  }
  return next(action);
};

let lastGCTime = 0;
const InterGCDelayMS = 30 * 1000; // 30 seconds
function _itIsTimeForGC() {
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

function _doGC(store) {
  const state = store.getState();

  // these should probably be a function imported from undoable.js, etc, as
  // they have overly intimiate knowledge of our reducers.
  const undoablePast = state["@@undoable/past"];
  const undoableFuture = state["@@undoable/future"];
  const undoableStack = undoablePast
    .concat(undoableFuture)
    .flatMap((snapshot) =>
      snapshot.filter((v) => v[0] === "annoMatrix").map((v) => v[1])
    );
  const currentAnnoMatrix = state.annoMatrix;

  /*
  We want to identify those matrixes currently "hot", ie, linked from the current annoMatrix,
  as our current gc algo is more aggressive with those not hot.
  */
  const allAnnoMatrices = new Map(
    undoableStack.map((m) => [m, { isHot: false }])
  );
  let am = currentAnnoMatrix;
  while (am) {
    allAnnoMatrices.set(am, { isHot: true });
    am = am.viewOf;
  }
  allAnnoMatrices.forEach((hints, annoMatrix) => annoMatrix._gc(hints));
}

export default annoMatrixGC;
