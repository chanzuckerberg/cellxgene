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
import { Reducer, AnyAction } from "redux";
import fromEntries from "../util/fromEntries";

export const pastKey = "@@undoable/past";
export const futureKey = "@@undoable/future";
export const filterStateKey = "@@undoable/filterState";
export const filterActionKey = "@@undoable/filterAction";
export const pendingKey = "@@undoable/pending";
const defaultHistoryLimit = -100;

export interface UndoableFilterState {
  [name: string]: unknown;
}

export interface UndoableConfig<FilterStateType extends UndoableFilterState> {
  debug?: boolean | number;
  historyLimit?: number;
  actionFilter?: ActionFilterFn<FilterStateType>;
}

export interface UndoableAction<FilterStateType extends UndoableFilterState> {
  [filterActionKey]: string;
  [filterStateKey]?: FilterStateType;
}

export type ActionFilterFn<FilterStateType extends UndoableFilterState> = (
  undoableState: UndoableState<FilterStateType>,
  action: AnyAction,
  filterState?: FilterStateType
) => UndoableAction<FilterStateType>;

export interface UndoableState<FilterStateType extends UndoableFilterState> {
  [pastKey]: [string, unknown][][];
  [futureKey]: [string, unknown][][];
  [pendingKey]: [string, unknown][] | null;
  [filterStateKey]: FilterStateType | undefined;
}

const Undoable = <FilterStateType extends UndoableFilterState>(
  reducer: Reducer,
  undoableKeys: string[],
  options: UndoableConfig<FilterStateType> = {}
): Reducer => {
  const debug = options?.debug ?? false;
  let { historyLimit } = options;
  if (!historyLimit) historyLimit = defaultHistoryLimit;
  if (historyLimit > 0) historyLimit = -historyLimit;
  const actionFilter: ActionFilterFn<FilterStateType> =
    options?.actionFilter ?? (() => ({ [filterActionKey]: "save" }));

  if (!Array.isArray(undoableKeys) || undoableKeys.length === 0)
    throw new Error("undoable keys array must be specified");
  const undoableKeysSet = new Set(undoableKeys);

  /*
  Undo the current to previous history
  */
  function undo(
    currentState: UndoableState<FilterStateType>
  ): UndoableState<FilterStateType> {
    const past = currentState[pastKey];
    const future = currentState[futureKey];
    if (past.length === 0) return currentState;
    const currentUndoableState = Object.entries(currentState).filter((kv) =>
      undoableKeysSet.has(kv[0])
    );
    const newPast = [...past];
    const newState = newPast.pop() || [];
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
  function redo(
    currentState: UndoableState<FilterStateType>
  ): UndoableState<FilterStateType> {
    const past = currentState[pastKey] || [];
    const future = currentState[futureKey] || [];
    if (future.length === 0) return currentState;
    const currentUndoableState = Object.entries(currentState).filter((kv) =>
      undoableKeysSet.has(kv[0])
    );
    const newFuture = [...future];
    const newState = newFuture.pop() || [];
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
  function clear(
    currentState: UndoableState<FilterStateType>
  ): UndoableState<FilterStateType> {
    return {
      ...currentState,
      [pastKey]: [],
      [futureKey]: [],
      [filterStateKey]: undefined,
      [pendingKey]: null,
    };
  }

  /*
  Reduce current action, with no history side-effects
  */
  function skip(
    currentState: UndoableState<FilterStateType>,
    action: AnyAction,
    filterState: UndoableFilterState
  ): UndoableState<FilterStateType> {
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
  function save(
    currentState: UndoableState<FilterStateType>,
    action: AnyAction,
    filterState: UndoableFilterState
  ): UndoableState<FilterStateType> {
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
  function stashPending(
    currentState: UndoableState<FilterStateType>
  ): UndoableState<FilterStateType> {
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
  function cancelPending(
    currentState: UndoableState<FilterStateType>
  ): UndoableState<FilterStateType> {
    return {
      ...currentState,
      [pendingKey]: null,
    };
  }

  /*
  Push pending state onto the history stack
  */
  function applyPending(
    currentState: UndoableState<FilterStateType>
  ): UndoableState<FilterStateType> {
    const past = currentState[pastKey];
    const pendingState = currentState[pendingKey];
    if (pendingState === null) return currentState;
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
    currentState: UndoableState<FilterStateType> = {
      [pastKey]: [],
      [futureKey]: [],
      [filterStateKey]: undefined,
      [pendingKey]: null,
    },
    action: AnyAction
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

function push<T = unknown>(arr: T[], val: T, limit?: number) {
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
