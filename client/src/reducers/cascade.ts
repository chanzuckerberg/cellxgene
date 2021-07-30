// eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any -- - FIXME: disabled temporarily on migrate to TS.
export default function cascadeReducers(arg: any) {
  /*
  Combine a set of cascading reducers into a single reducer.  Cascading
  reducers are reducers which may rely on state computed by another reducer.
  Therefore, they:
  - must be composed in a particular order (currently, this is a simple
    linear list of reducers, run in list order)
  - must have access to partially updated "next state" so they can further
    derive state.

  Parameter is one of:
  - a Map object
  - an array of tuples, [ [key1, reducer1], [key2, reducer2], ... ]
  Ie,  cascadeReducers([ ["a", reduceA], ["b", reduceB] ])

  Each reducer will be called with the signature:
      (prevState, action, sharedNextState, sharedPrevState) => newState

  cascadeReducers will build a composite newState object, much
  like combinedReducers.  Additional semantics:
  - reducers guaranteed to be called in order
  - each reducer will receive shared objects
  */
  const reducers = arg instanceof Map ? arg : new Map(arg);
  const reducerKeys = [...reducers.keys()];
  // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
  return (prevState: any, action: any) => {
    const nextState = {};
    let stateChange = false;
    for (let i = 0, l = reducerKeys.length; i < l; i += 1) {
      const key = reducerKeys[i];
      const reducer = reducers.get(key);
      const prevStateForKey = prevState ? prevState[key] : undefined;
      const nextStateForKey = reducer(
        prevStateForKey,
        action,
        nextState,
        prevState
      );
      // @ts-expect-error ts-migrate(7053) FIXME: Element implicitly has an 'any' type because expre... Remove this comment to see the full error message
      nextState[key] = nextStateForKey;
      stateChange = stateChange || nextStateForKey !== prevStateForKey;
    }
    return stateChange ? nextState : prevState;
  };
}
