/*
Garbage collection / cache management support

Middleware that knows how to garbage collect undo state.
*/

const annoMatrixGC = (store) => (next) => (action) => {
  if (_itIsTimeForGC()) {
    _doGC(store);
  }
  return next(action);
};

let lastGCTime = 0;
const ThirtySeconds = 30 * 1000;
function _itIsTimeForGC() {
  /*
  we don't want to run GC on every dispatch, so throttle it a bit.
  */
  const now = Date.now();
  if (now - lastGCTime > ThirtySeconds) {
    lastGCTime = now;
    return true;
  }
  return false;
}

function _doGC(store) {
  console.log("...running GC middleware");

  const state = store.getState();

  // these should probably be a function imported from undoable.js, etc, as
  // they have overly intimiate knowledge of our reducers.
  const undoablePast = state["@@undoable/past"];
  const undoableFuture = state["@@undoable/future"];
  const undoableStack = undoablePast.concat(undoableFuture).flatMap((snapshot) =>
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
