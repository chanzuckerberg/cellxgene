import StateMachine from "../util/statemachine";
import createFsmTransitions from "./undoableFsm";

const actionKey = "@@undoable/filterAction";
const stateKey = "@@undoable/filterState";

/*
these actions will not affect history
*/
const skipOnActions = new Set([
  "url changed",
  "interface reset started",
  "initial data load start",
  "configuration load complete",
  "increment graph render counter",
  "window resize",
  "user reset start",
  "reset colorscale",

  "graph brush change",
  "continuous metadata histogram brush",

  "expression load start",
  "expression load success",
  "expression load error",

  "request user defined gene started",
  "request user defined gene success",
  "clear all user defined genes",

  "get single gene expression for coloring started",
  "get single gene expression for coloring error"
]);

/*
identical, repeated occurances of these action types will be debounced.
Entire action must be identical (all keys).
*/
const debounceOnActions = new Set([
  "color by categorical metadata",
  "color by continuous metadata",
  "color by expression"
]);

/*
history will be cleared when these actions occur
*/
const clearOnActions = new Set([
  "initial data load complete (universe exists)",
  "reset World to eq Universe",
  "initial data load error",
  "user reset end"
]);

/*
An immediate history save will be done for these
*/
const saveOnActions = new Set([
  "categorical metadata filter select",
  "categorical metadata filter deselect",
  "categorical metadata filter all of these",
  "categorical metadata none of these",

  "color by categorical metadata",
  "color by continuous metadata",
  "color by expression",

  "set scatterplot x",
  "set scatterplot y",

  "store current cell selection as differential set 1",
  "store current cell selection as differential set 2",

  "set World to current selection"
]);

/**
StateMachine - processing complex action handling - see FSM graph for
actual structure, in undoableFsm.js
**/

/*
Default FSM actions.  Used to side-effect transitions in the graph.
See graph definition for the transitions that use each.

Signature:  (fsm, transition, reducerState, reducerAction) => undoableAction
*/
const stashPending = fsm => ({
  [actionKey]: "stashPending",
  [stateKey]: { fsm }
});
const cancelPending = () => ({
  [actionKey]: "cancelPending",
  [stateKey]: { fsm: null }
});
const applyPending = () => ({
  [actionKey]: "applyPending",
  [stateKey]: { fsm: null }
});
const skip = fsm => ({ [actionKey]: "skip", [stateKey]: { fsm } });
const clear = () => ({ [actionKey]: "clear", [stateKey]: { fsm: null } });
const save = fsm => ({ [actionKey]: "save", [stateKey]: { fsm } });

/*
Error handler for state transitions that are unexpected.  Called by
StateMachine when it doesn't know what to do.

Signature:  (fsm, event, from) => undoableAction
*/
const onFsmError = (fsm, name, from) => {
  console.error("FSM error - unexpected history state", fsm, name, from);
  // In production, try to recover gracefully if we have unexpected state
  return clear(fsm);
};

/*
Definition of the transition graph mapping action types to history side effects.
*/
const fsmTransitions = createFsmTransitions(
  stashPending,
  cancelPending,
  applyPending,
  skip,
  clear,
  save
);
/* State machine we clone whenever we need to run it */
const seedFsm = new StateMachine("init", fsmTransitions, onFsmError);

/*
See undoable.js for description action filter interface description.

Basic approach:
  * trivial handlers for skip, clear & save cases to keep config simple.
  * only implement complex state machines where absolutely required (eg,
    multi-event seleciton and the like)
*/
const actionFilter = debug => (state, action, prevFilterState) => {
  const actionType = action.type;
  const filterState = {
    ...prevFilterState,
    prevAction: action
  };
  if (skipOnActions.has(actionType)) {
    return { [actionKey]: "skip", [stateKey]: filterState };
  }
  if (
    debounceOnActions.has(actionType) &&
    shallowObjectEq(action, prevFilterState.prevAction)
  ) {
    return { [actionKey]: "skip", [stateKey]: filterState };
  }
  if (clearOnActions.has(actionType)) {
    return { [actionKey]: "clear", [stateKey]: filterState };
  }
  if (saveOnActions.has(actionType)) {
    return { [actionKey]: "save", [stateKey]: filterState };
  }

  /*
    Else, something more complex OR unknown to us....
    */
  if (seedFsm.events.has(actionType)) {
    let { fsm } = filterState;
    if (!fsm) {
      /* no active FSM, so create one in init state */
      fsm = seedFsm.clone("init");
    }
    return fsm.next(action.type, { state, action });
  }

  /* else, we have no idea what this is - skip it */
  if (debug) console.log("**** ACTION FILTER EVENT HANDLER MISS", actionType);
  return { [actionKey]: "skip", [stateKey]: filterState };
};

/*
return true if objA and objB are ===, OR if:
  - are both objects and not null
  - have same own properties
  - all values are strict equal (===)
*/
function shallowObjectEq(objA, objB) {
  if (objA === objB) return true;
  if (!objA || !objB) return false;
  if (!shallowArrayEq(Object.keys(objA), Object.keys(objB))) return false;
  if (!shallowArrayEq(Object.values(objA), Object.values(objB))) return false;
  return true;
}

/*
return true if arrA and arrB contain the same strict-equal values,
in the same order.
*/
function shallowArrayEq(arrA, arrB) {
  if (arrA.length !== arrB.length) return false;
  for (let i = 0, l = arrA.length; i < l; i += 1) {
    if (arrA[i] !== arrB[i]) return false;
  }
  return true;
}

/* configuration for the undoable meta reducer */
const debug = false; // set truish for undoble debugging
const undoableConfig = {
  debug,
  historyLimit: 50, // maximum history size
  actionFilter: actionFilter(debug)
};

/*
this code is strictly for sanity checking configuration, and is only
enabled when we are debugging the undoable configuration (ie, debug === true).
*/
if (debug) {
  /*
  Confirm no intersection between the various trivial rejection action filters
  */
  if (
    new Set([...skipOnActions].filter(x => clearOnActions.has(x))).size > 0 ||
    new Set([...skipOnActions].filter(x => saveOnActions.has(x))).size > 0 ||
    new Set([...clearOnActions].filter(x => saveOnActions.has(x))).size > 0
  ) {
    console.error(
      "Undoable misconfiguration - action filters have redundant events"
    );
  }

  /*
  Confirm that no FSM events are blocked by a trivial rejection filter.
  If this occurs, the FSM can't ever see the events needed to process
  state transitions.
  */
  const trivialFilters = new Set([
    ...skipOnActions,
    ...clearOnActions,
    ...saveOnActions
  ]);
  const trivialOverlapWithFsm = new Set(
    [...trivialFilters].filter(x => seedFsm.events.has(x))
  );
  if (trivialOverlapWithFsm.size > 0) {
    console.error(
      "Undoable misconfiguration - trivival action filter blocking FSM filter",
      [...trivialOverlapWithFsm]
    );
  }
}

export default undoableConfig;
