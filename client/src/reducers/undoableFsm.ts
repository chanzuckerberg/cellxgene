/*
State transition graph for complex action/history interactions.

Assumed configuration from undoableConfig:
    * By convention, "init" is used as the start state for all, and "done"
      as the final state.
    * Unexpected states will result in an error, plus a clear and cancelPending
      side-effect.

TODO: is is possible there is a more concise format for this, as it is
a fairly repetitive pattern.

These events are largely one of two types:
a) async operations or multi-event options that should only be committed
upon some success criteria, otherwise cancelled.

b) compound actions that should be collapsed into a single history change.

*/

// eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types --- FIXME: disabled temporarily on migrate to TS.
const createFsmTransitions = (
  // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any -- - FIXME: disabled temporarily on migrate to TS.
  stashPending: any,
  // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any -- - FIXME: disabled temporarily on migrate to TS.
  cancelPending: any,
  // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any -- - FIXME: disabled temporarily on migrate to TS.
  applyPending: any,
  // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any -- - FIXME: disabled temporarily on migrate to TS.
  skip: any,
  // @ts-expect-error ts-migrate(6133) FIXME: 'clear' is declared but its value is never read.
  // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any -- - FIXME: disabled temporarily on migrate to TS.
  clear: any,
  // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any -- - FIXME: disabled temporarily on migrate to TS.
  save: any
) => {
  return [
    /* graph selection brushing */
    {
      event: "graph brush start",
      from: "init",
      to: "graph brush in progress",
      action: stashPending,
    },
    {
      event: "graph brush cancel",
      from: "graph brush in progress",
      to: "done",
      action: applyPending,
    },
    {
      event: "graph brush deselect",
      from: "graph brush in progress",
      to: "done",
      /* if current selection is all, cancelPending.  Else, applyPending */
      // @ts-expect-error ts-migrate(6133) FIXME: 'fsm' is declared but its value is never read.
      // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
      action: (fsm: any, transition: any, data: any) =>
        data.state.graphSelection.selection.mode === "all"
          ? cancelPending()
          : applyPending(),
    },
    {
      event: "graph brush end",
      from: "graph brush in progress",
      to: "done",
      action: applyPending,
    },

    /* graph selection lasso */
    {
      event: "graph lasso start",
      from: "init",
      to: "graph lasso in progress",
      action: stashPending,
    },
    {
      event: "graph lasso cancel",
      from: "graph lasso in progress",
      to: "done",
      action: applyPending,
    },
    {
      event: "graph lasso deselect",
      from: "graph lasso in progress",
      to: "done",
      /* if current selection is all, cancelPending.  Else, applyPending */
      // @ts-expect-error ts-migrate(6133) FIXME: 'fsm' is declared but its value is never read.
      // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
      action: (fsm: any, transition: any, data: any) =>
        data.state.graphSelection.selection.mode === "all"
          ? cancelPending()
          : applyPending(),
    },
    {
      event: "graph lasso end",
      from: "graph lasso in progress",
      to: "done",
      action: applyPending,
    },

    /* Continuous metadata histogram brush selection */
    {
      event: "continuous metadata histogram start",
      from: "init",
      to: "continuous histo select in progress",
      action: stashPending,
    },
    {
      event: "continuous metadata histogram cancel",
      from: "continuous histo select in progress",
      to: "done",
      action: cancelPending,
    },
    {
      event: "continuous metadata histogram cancel",
      from: "init",
      to: "done",
      action: save,
    },
    {
      event: "continuous metadata histogram end",
      from: "continuous histo select in progress",
      to: "done",
      action: applyPending,
    },

    /* Single gene request by user */
    {
      event: "single user defined gene start",
      from: "init",
      to: "single user gene request in progress",
      action: stashPending,
    },
    {
      event: "request user defined gene error",
      from: "single user gene request in progress",
      to: "single user gene error in progress",
      action: skip,
    },
    {
      event: "single user defined gene error",
      from: "single user gene error in progress",
      to: "done",
      action: cancelPending,
    },
    {
      event: "single user defined gene complete",
      from: "single user gene request in progress",
      to: "done",
      action: applyPending,
    },

    /* Bulk gene request by user */
    {
      event: "bulk user defined gene start",
      from: "init",
      to: "bulk user gene request in progress",
      action: stashPending,
    },
    {
      event: "request user defined gene error",
      from: "bulk user gene request in progress",
      to: "bulk user gene request error in progress",
      action: skip,
    },
    {
      event: "bulk user defined gene error",
      from: "bulk user gene request error in progress",
      to: "done",
      action: cancelPending,
    },
    {
      event: "bulk user defined gene complete",
      from: "bulk user gene request in progress",
      to: "done",
      action: applyPending,
    },

    /* Compute Differential Expression button user action */
    {
      event: "request differential expression started",
      from: "init",
      to: "diffexp in progress",
      action: stashPending,
    },
    {
      event: "request user defined gene error",
      from: "diffexp in progress",
      to: "done",
      action: cancelPending,
    },
    {
      event: "request differential expression success",
      from: "diffexp in progress",
      to: "done",
      action: applyPending,
    },

    /* clear scatter plot button (eg, on scatterplot view) */
    {
      event: "clear scatterplot",
      from: "init",
      to: "done",
      action: save,
    },
  ];
};

export default createFsmTransitions;
