/*
A redo/undo meta reducer for Redux.  Designed to work well with the cascadeReducer().

Requires three parameters:
* reducer - a reducer, which MUST return an object as state.
* undoableKeys - an array of object keys (strings).  If any of these keys
  are in the object/state returned by the reducer, they will be treated as
  state to be made "undoable".
* options - an optional object, which may contain the following parameters:
    * historyLimit: max number of historical states to remember (aka max undo depth)
    * actionFilter: filter function, (state, action, filterState) => value.
      See below for details.
    * debug: if truish, will print helpful log messages about history manipulation

This meta reducer accepts three actions types:
* @@undoable/undo - move back in history
* @@undoable/redo - move forward in history
* @@undoable/clear - clear history

---

Action filter - controls the undoable reducer side-effects.  If not
specified, the action filter defaults to "save", ie, pushes a redo
point upon each action.

The action filter callback has access to the current action, the entire
undoable reducer state, and any state it wants to manage ("filterState").
This filter state will be passed to each action filter call, and any
value returned (via @@undoable/filterState field described below) will be
MERGED into the current filter state.

An object must be returned (the "undoable action"), indicating desired
history state processing.  The undoable action object contents, by key:

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
    MERGED into the current filter state. The value and semantics of any
    filter state are entirely at the discretion of the action filter.

*/
import fromEntries from "../util/fromEntries";

const historyKeyPrefix = "@@undoable/";
const pastKey = `${historyKeyPrefix}past`;
const futureKey = `${historyKeyPrefix}future`;
const filterStateKey = `${historyKeyPrefix}filterState`;
const filterActionKey = `${historyKeyPrefix}filterAction`;
const pendingKey = `${historyKeyPrefix}pending`;
const defaultHistoryLimit = -100;

// eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any -- - FIXME: disabled temporarily on migrate to TS.
const Undoable = (reducer: any, undoableKeys: any, options = {}) => {
  // @ts-expect-error ts-migrate(2339) FIXME: Property 'debug' does not exist on type '{}'.
  const { debug } = options;
  // @ts-expect-error ts-migrate(2339) FIXME: Property 'historyLimit' does not exist on type '{}... Remove this comment to see the full error message
  let { historyLimit } = options;
  if (!historyLimit) historyLimit = defaultHistoryLimit;
  if (historyLimit > 0) historyLimit = -historyLimit;
  const actionFilter =
    // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
    (options as any).actionFilter || (() => ({ [filterActionKey]: "save" }));

  if (!Array.isArray(undoableKeys) || undoableKeys.length === 0)
    throw new Error("undoable keys array must be specified");
  const undoableKeysSet = new Set(undoableKeys);

  /*
  Undo the current to previous history
  */
  // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
  function undo(currentState: any) {
    const past = currentState[pastKey];
    const future = currentState[futureKey];
    if (past.length === 0) return currentState;
    const currentUndoableState = Object.entries(currentState).filter((kv) =>
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
      [pendingKey]: null,
    };
    return nextState;
  }

  /*
  Replay future, previously undone.
  */
  // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
  function redo(currentState: any) {
    const past = currentState[pastKey] || [];
    const future = currentState[futureKey] || [];
    if (future.length === 0) return currentState;
    const currentUndoableState = Object.entries(currentState).filter((kv) =>
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
      [pendingKey]: null,
    };
    return nextState;
  }

  /*
  Clear the history state.  No side-effects on current state.
  */
  // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
  function clear(currentState: any) {
    return {
      ...currentState,
      [pastKey]: [],
      [futureKey]: [],
      [filterStateKey]: {},
      [pendingKey]: null,
    };
  }

  /*
  Reduce current action, with no history side-effects
  */
  // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
  function skip(currentState: any, action: any, filterState: any) {
    const past = currentState[pastKey] || [];
    const future = currentState[futureKey] || [];
    const pending = currentState[pendingKey];
    const res = reducer(currentState, action);
    return {
      ...res,
      [pastKey]: past,
      [futureKey]: future,
      [filterStateKey]: filterState,
      [pendingKey]: pending,
    };
  }

  /*
  Save current state in the history, then reduce action.
  */
  // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
  function save(currentState: any, action: any, filterState: any) {
    const past = currentState[pastKey] || [];
    const currentUndoableState = Object.entries(currentState).filter((kv) =>
      undoableKeysSet.has(kv[0])
    );
    const res = reducer(currentState, action);
    const newPast = push(past, currentUndoableState, historyLimit);
    const nextState = {
      ...res,
      [pastKey]: newPast,
      [futureKey]: [],
      [filterStateKey]: filterState,
      [pendingKey]: null,
    };
    return nextState;
  }

  /*
  Save current state as pending history change.  No other side effects.
  */
  // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
  function stashPending(currentState: any) {
    const currentUndoableState = Object.entries(currentState).filter((kv) =>
      undoableKeysSet.has(kv[0])
    );
    return {
      ...currentState,
      [pendingKey]: currentUndoableState,
    };
  }

  /*
  Cancel pending history state change.  No other side effects.
  */
  // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
  function cancelPending(currentState: any) {
    return {
      ...currentState,
      [pendingKey]: null,
    };
  }

  /*
  Push pending state onto the history stack
  */
  // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
  function applyPending(currentState: any) {
    const past = currentState[pastKey] || [];
    const pendingState = currentState[pendingKey];
    const newPast = push(past, pendingState, historyLimit);
    const nextState = {
      ...currentState,
      [pastKey]: newPast,
      [futureKey]: [],
      [pendingKey]: null,
    };
    return nextState;
  }

  return (
    currentState = {
      [pastKey]: [],
      [futureKey]: [],
      [filterStateKey]: {},
      [pendingKey]: null,
    },
    // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
    action: any
  ) => {
    if (debug > 1) console.log("---- ACTION", action.type);
    const aType = action.type;
    switch (aType) {
      case "@@undoable/undo": {
        return undo(currentState);
      }

      case "@@undoable/redo": {
        return redo(currentState);
      }

      case "@@undoable/clear": {
        return clear(currentState);
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
          [filterStateKey]: filterStateUpdate,
        } = actionFilterResp;
        const nextFilterState = { ...currentFilterState, ...filterStateUpdate };

        switch (filterAction) {
          case "clear":
            if (debug) console.log("---- CLEAR HISTO", action.type);
            return clear(skip(currentState, action, nextFilterState));

          case "save":
            if (debug) console.log("---- SAVE HISTO", action.type);
            return save(currentState, action, nextFilterState);

          case "stashPending":
            if (debug) console.log("---- STASH PENDING", action.type);
            return skip(stashPending(currentState), action, nextFilterState);

          case "cancelPending":
            if (debug) console.log("---- CANCEL PENDING", action.type);
            return skip(cancelPending(currentState), action, nextFilterState);

          case "applyPending":
            if (debug) console.log("---- APPLY PENDING", action.type);
            return skip(applyPending(currentState), action, nextFilterState);

          case "skip":
          default:
            return skip(currentState, action, nextFilterState);
        }
      }
    }
  };
};

// eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
function push(arr: any, val: any, limit = undefined) {
  /*
  functional array push, with a max length limit to the new array.
  Like Array.push, except it returns new array and discards as needed
  to enforce the length limit.
  */
  const narr = arr.slice(limit);
  narr.push(val);
  return narr;
}

export default Undoable;
