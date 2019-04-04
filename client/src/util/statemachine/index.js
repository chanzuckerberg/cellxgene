/*
Very simple FSM for use in reducer, etc.
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
    this.events = new Set(transitions.map(tsn => tsn.event));

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

  next(event) {
    const { graph, state } = this;
    const tsnMap = graph.get(event);
    if (!tsnMap) return this.onError(this, event, state, undefined);

    const transition = tsnMap.get(state);
    if (!transition) return this.onError(this, event, state, undefined);

    this.state = transition.to;
    return transition.action ? transition.action(this, transition) : undefined;
  }
}
