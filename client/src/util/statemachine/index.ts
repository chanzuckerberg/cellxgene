/*
Very simple FSM for use in reducer, etc.

To create a state machine:
  new StateMachine(initialState, transitions, onErrorCallback) -> statemachine

Where:
  * initialState - a caller-specified value that represents the initial state of
    the FSM.
  * transitions - an array of objects, representing FSM transitions (graph edges),
    having the form:
      {
        to: state_name_transitioning_to,
        from: state_name_transitioning_from,
        event: value_that_will_cause_transition,
        action: optional_callback_upon_transition
      }
    The transition will be provided to the action callback, so other data
    may be stored in the transition object for use by the action callback.
  * onErrorCallback - a callback function called if the FSM receives an event
    for which it has no defined transition.


Interface:
  * states - property containing the state names. A Set(), contianing the
    union of to: and from: values.
  * events - property containing all of the accepted event values.  Set().
  * graph - a Map of Maps, organized as  graph[eventValue][fromStateValue]
  * clone() - clone the entire statemachine.
  * next(eventValue) - drive the FSM to the next state.   If the event
    matches a transition with a defined action, the action callback is
    called, and the action return value is returned by next().  If no
    transition is defined, onErrorCallback is called.

Example:

  const transitions = [
    { from: "A", to: "B", event: "yo", action: () => 42 }
  ];
  const fsm = new StateMachine("A", transitions, () => { throw new Error("oops") });
  fsm.next("yo");   // returns 42


*/
export default class StateMachine {
  constructor(initState, transitions, onError) {
    this.onError = onError || (() => undefined);
    this.state = initState;

    // all states
    this.states = new Set(
      transitions.reduce((names, tsn) => {
        names.push(tsn.from);
        names.push(tsn.to);
        return names;
      }, [])
    );

    // all transition names (aka events)
    this.events = new Set(transitions.map((tsn) => tsn.event));

    // the transition graph.
    // graph[event][from] -> transition
    this.graph = transitions.reduce((graph, tsn) => {
      const { event, from } = tsn;
      if (!graph.has(event)) graph.set(event, new Map());
      const tsnMap = graph.get(event);
      tsnMap.set(from, tsn);
      return graph;
    }, new Map());
  }

  clone(initState) {
    const fsm = new StateMachine(initState, []);
    fsm.onError = this.onError;
    fsm.states = this.states;
    fsm.events = this.events;
    fsm.graph = this.graph;
    return fsm;
  }

  next(event, data) {
    const { graph, state } = this;
    const tsnMap = graph.get(event);
    if (!tsnMap) return this.onError(this, event, state, undefined);

    const transition = tsnMap.get(state);
    if (!transition) return this.onError(this, event, state, undefined);

    this.state = transition.to;
    return transition.action
      ? transition.action(this, transition, data)
      : undefined;
  }
}
