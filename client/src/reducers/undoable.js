/*
A redo/undo meta reducer for Redux.   Designed to work well with the cascadeReducer().

Requires three parameters:
* reducer - a reducer, which MUST return an object as state.
* undoableKeys - an array of object keys (strings).   If any of these keys
  are in the object/state returned by the reducer, they will be treated as
  state to be made "undoable".
* options - an optional object, which may contain the following parameters:
    * historyLimit: max number of historical states to remember (aka max undo depth)
    * skipActionFilter: filter function, (state, action) => bool.  If it returns
      truthy, the current state will not be pushed onto the history stack.
    * clearOnActionFilter: filter function, (state, action) => bool.  If it returns
      truthy, the history state will be cleared as part of handling this action.

skipActionFilter has precedence over clearOnActionFilter.

This meta reducer accepts three actions types:
* @@undoable/undo - move back in history
* @@undoable/redo - move forward in history
* @@undoable/clear - clear history
*/

const historyKeyPrefix = "@@undoable/";
const pastKey = `${historyKeyPrefix}past`;
const futureKey = `${historyKeyPrefix}future`;
const defaultHistoryLimit = -100;

const Undoable = (reducer, undoableKeys, options = {}) => {
  let { historyLimit } = options;
  if (!historyLimit) historyLimit = defaultHistoryLimit;
  if (historyLimit > 0) historyLimit = -historyLimit;
  const skipActionFilter = options.skipActionFilter || (() => false);
  const clearOnActionFilter = options.clearOnActionFilter || (() => false);

  if (!Array.isArray(undoableKeys) || undoableKeys.length === 0)
    throw new Error("undoable keys array must be specified");
  const undoableKeysSet = new Set(undoableKeys);

  function undo(currentState) {
    const past = currentState[pastKey];
    const future = currentState[futureKey];
    if (past.length === 0) return currentState;
    const currentUndoableState = Object.entries(currentState).filter(kv =>
      undoableKeysSet.has(kv[0])
    );
    const newPast = [...past];
    const newState = newPast.pop();
    const newFuture = push(future, currentUndoableState);
    const nextState = {
      ...currentState,
      ...fromEntries(newState),
      [pastKey]: newPast,
      [futureKey]: newFuture
    };
    return nextState;
  }

  function redo(currentState) {
    const past = currentState[pastKey] || [];
    const future = currentState[futureKey] || [];
    if (future.length === 0) return currentState;
    const currentUndoableState = Object.entries(currentState).filter(kv =>
      undoableKeysSet.has(kv[0])
    );
    const newFuture = [...future];
    const newState = newFuture.pop();
    const newPast = push(past, currentUndoableState);
    const nextState = {
      ...currentState,
      ...fromEntries(newState),
      [pastKey]: newPast,
      [futureKey]: newFuture
    };
    return nextState;
  }

  function clear(currentState) {
    return {
      ...currentState,
      [pastKey]: [],
      [futureKey]: []
    };
  }

  function skip(currentState, action) {
    const past = currentState[pastKey] || [];
    const res = reducer(currentState, action);
    return {
      ...res,
      [pastKey]: past,
      [futureKey]: []
    };
  }

  function save(currentState, action) {
    const past = currentState[pastKey] || [];
    const currentUndoableState = Object.entries(currentState).filter(kv =>
      undoableKeysSet.has(kv[0])
    );
    const res = reducer(currentState, action);
    const newPast = push(past, currentUndoableState, historyLimit);
    const nextState = {
      ...res,
      [pastKey]: newPast,
      [futureKey]: []
    };
    return nextState;
  }

  return (
    currentState = {
      [pastKey]: [],
      [futureKey]: []
    },
    action
  ) => {
    const aType = action.type;
    switch (aType) {
      case "@@undoable/undo": {
        return undo(currentState, action);
      }
      case "@@undoable/redo": {
        return redo(currentState, action);
      }
      case "@@undoable/clear": {
        return clear(currentState, action);
      }
      default: {
        if (skipActionFilter(currentState, action)) {
          return skip(currentState, action);
        }
        if (clearOnActionFilter(currentState, action)) {
          return clear(skip(currentState, action));
        }
        return save(currentState, action);
      }
    }
  };
};

function push(arr, val, limit = undefined) {
  /*
  functional array push, with a max length limit to the new array.
  Like Array.push, except it returns new array and discards as needed
  to enforce the length limit.
  */
  const narr = arr.slice(limit);
  narr.push(val);
  return narr;
}

function fromEntries(arr) {
  /*
  Similar to Object.fromEntries, but only handles array.
  This could be replaced with the standard fucnction once it
  is widely available.   As of 3/20/2019, it has not yet
  been released in the Chrome stable channel.
  */
  const obj = {};
  for (let i = 0, l = arr.length; i < l; i += 1) {
    obj[arr[i][0]] = arr[i][1];
  }
  return obj;
}

export default Undoable;
