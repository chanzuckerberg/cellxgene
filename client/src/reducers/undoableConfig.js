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

  "reset colorscale",
  "color by categorical metadata",
  "color by continuous metadata",
  "color by expression",

  "set scatterplot x",
  "set scatterplot y",
  "clear scatterplot",

  "clear differential expression",
  "store current cell selection as differential set 1",
  "store current cell selection as differential set 2",

  "set World to current selection"
]);

/*
StateMachine processed / complex action handling - see FSM graph for
actual structure.
*/

/* action args:   fsm, transition, reducerState, reducerAction */
const stashPending = fsm => ({ [actionKey]: "stashPending", [stateKey]: fsm });
const cancelPending = () => ({
  [actionKey]: "cancelPending",
  [stateKey]: null
});
const applyPending = () => ({ [actionKey]: "applyPending", [stateKey]: null });
const skip = fsm => ({ [actionKey]: "skip", [stateKey]: fsm });
const clear = () => ({ [actionKey]: "clear", [stateKey]: null });
const save = fsm => ({ [actionKey]: "save", [stateKey]: fsm });

/* onFsmError args:  fsm, name, from */
const onFsmError = (fsm, name, from) => {
  console.error("FSM error - unexpected history state", fsm, name, from);
  // In production, try to recover gracefully if we have unexpected state
  return clear(fsm);
};

/*
Definition of the transition graph mapping action types to history side effects
*/
const fsmTransitions = createFsmTransitions(
  stashPending,
  cancelPending,
  applyPending,
  skip,
  clear,
  save
);
const seedFsm = new StateMachine("init", fsmTransitions, onFsmError);

/*
See undoable.js for description action filter interface description

Basic approach:
  * trivial handlers for skip, clear & save cases to keep config simple.
  * only implement complex state machines where absolutely required (eg,
    multi-event seleciton and the like)
*/

function actionFilter(state, action, fsm) {
  if (skipOnActions.has(action.type)) {
    return skip(fsm);
  }
  if (clearOnActions.has(action.type)) {
    return clear(fsm);
  }
  if (saveOnActions.has(action.type)) {
    return save(fsm);
  }

  /*
  Else, something more complex OR unknown to us....
  */
  if (seedFsm.events.has(action.type)) {
    if (!fsm) {
      /* no active FSM, so create one in init state */
      fsm = seedFsm.clone("init");
    }
    return fsm.next(action.type, state, action);
  }

  /* else, we have no idea what this is - skip it */
  console.log("**** EVENT MISS", action.type);
  return skip(fsm);
}

/* configuration for the undoable meta reducer */
const undoableConfig = {
  historyLimit: 50, // maximum history size
  actionFilter
};

export default undoableConfig;
