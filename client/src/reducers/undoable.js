/*
A redo/undo meta reducer for Redux.   Designed to work well with the cascadeReducer().

Requires three parameters:
* reducer - a reducer, which MUST return an object as state.
* undoableKeys - an array of object keys (strings).   If any of these keys
  are in the object/state returned by the reducer, they will be treated as
  state to be made "undoable".
* options - an optional object, which may contain the following parameters:
    * historyLimit: max number of historical states to remember (aka max undo depth)
    * actionFilter: filter function, (state, action) => value.  See below for
      details

This meta reducer accepts three actions types:
* @@undoable/undo - move back in history
* @@undoable/redo - move forward in history
* @@undoable/clear - clear history


Action filter - controls the undoable reducer side-effects.
This will have access to both the current action, and the entire
undoable reducer state.

An object must be returned, indicating desired history state processing.
The object keys:

  @@undoable/filterAction: required.   Can be one of:
      "skip" - reduce the current action, but no other side effects.
        Same as returning false.
      "clear" - reduce the current action, and clear history state.
      "save" - push the previous state onto the history stack (ie,
        before reducing the action)
      "stashPending" - reduce action, save state as pending.  Does not
        not commit it to history.  Along with cancelPending and applyPending,
        can be used to delay commit of history (eg, for multi-action
        groupings, asynch operations, etc).
      "cancelPending" - reduce action, cancel any pending state save.
      "applyPending" - commit any pending state to the history stack,
        then reduce action.

  @@undoable/filterState: optional.  If this value is set, it will be
    saved into the state field of the same name. The value and its
    semantics are entirely defined by the Action Filter function.

*/

const historyKeyPrefix = "@@undoable/";
const pastKey = `${historyKeyPrefix}past`;
const futureKey = `${historyKeyPrefix}future`;
const filterStateKey = `${historyKeyPrefix}filterState`;
const filterActionKey = `${historyKeyPrefix}filterAction`;
const pendingKey = `${historyKeyPrefix}pending`;
const defaultHistoryLimit = -100;

const Undoable = (reducer, undoableKeys, options = {}) => {
  let { historyLimit } = options;
  if (!historyLimit) historyLimit = defaultHistoryLimit;
  if (historyLimit > 0) historyLimit = -historyLimit;
  const actionFilter =
    options.actionFilter || (() => ({ [filterActionKey]: "save" }));

  if (!Array.isArray(undoableKeys) || undoableKeys.length === 0)
    throw new Error("undoable keys array must be specified");
  const undoableKeysSet = new Set(undoableKeys);

  /*
  Undo the current to previous history
  */
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
      [futureKey]: newFuture,
      [pendingKey]: null
    };
    return nextState;
  }

  /*
  Replay future, previously undone.
  */
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
      [futureKey]: newFuture,
      [pendingKey]: null
    };
    return nextState;
  }

  /*
  Clear the history state.  No side-effects on current state.
  */
  function clear(currentState) {
    return {
      ...currentState,
      [pastKey]: [],
      [futureKey]: [],
      [filterStateKey]: null,
      [pendingKey]: null
    };
  }

  /*
  Reduce current action, with no history side-effects
  */
  function skip(currentState, action, filterState) {
    const past = currentState[pastKey] || [];
    const pending = currentState[pendingKey];
    const res = reducer(currentState, action);
    return {
      ...res,
      [pastKey]: past,
      [futureKey]: [],
      [filterStateKey]: filterState,
      [pendingKey]: pending
    };
  }

  /*
  Save current state in the history, then reduce action.
  */
  function save(currentState, action, filterState) {
    const past = currentState[pastKey] || [];
    const currentUndoableState = Object.entries(currentState).filter(kv =>
      undoableKeysSet.has(kv[0])
    );
    const res = reducer(currentState, action);
    const newPast = push(past, currentUndoableState, historyLimit);
    const nextState = {
      ...res,
      [pastKey]: newPast,
      [futureKey]: [],
      [filterStateKey]: filterState,
      [pendingKey]: null
    };
    return nextState;
  }

  /*
  Save current state as pending history change.  No other side effects.
  */
  function stashPending(currentState) {
    const currentUndoableState = Object.entries(currentState).filter(kv =>
      undoableKeysSet.has(kv[0])
    );
    return {
      ...currentState,
      [pendingKey]: currentUndoableState
    };
  }

  /*
  Cancel pending history state change.  No other side effects.
  */
  function cancelPending(currentState) {
    return {
      ...currentState,
      [pendingKey]: null
    };
  }

  /*
  Push pending state onto the history stack
  */
  function applyPending(currentState) {
    const past = currentState[pastKey] || [];
    const pendingState = currentState[pendingKey];
    const newPast = push(past, pendingState, historyLimit);
    const nextState = {
      ...currentState,
      [pastKey]: newPast,
      [futureKey]: [],
      [pendingKey]: null
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
        const currentFilterState = currentState[filterStateKey];
        const actionFilterResp = actionFilter(
          currentState,
          action,
          currentFilterState
        );
        const {
          [filterActionKey]: filterAction,
          [filterStateKey]: nextFilterState
        } = actionFilterResp;

        switch (filterAction) {
          case "clear":
            // console.log("---- CLEAR HISTO", action.type);
            return clear(skip(currentState, action, nextFilterState));

          case "save":
            // console.log("---- SAVE HISTO", action.type);
            return save(currentState, action, nextFilterState);

          case "stashPending":
            // console.log("---- STASH", action.type);
            return skip(stashPending(currentState), action, nextFilterState);

          case "cancelPending":
            // console.log("---- CANCEL PENDING HISTO", action.type);
            return skip(cancelPending(currentState), action, nextFilterState);

          case "applyPending":
            // console.log("---- APPLY PENDING HISTO", action.type);
            return skip(applyPending(currentState), action, nextFilterState);

          case "skip":
          default:
            return skip(currentState, action, nextFilterState);
        }
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
