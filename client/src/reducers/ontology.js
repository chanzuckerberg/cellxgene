// jshint esversion: 6
const Ontology = (
  state = {
    enabled: false, // are ontology terms enabled?
    terms: null, // an array of term names, eg, ['cell', 'lung cell', ...]
    loading: true
  },
  action
) => {
  switch (action.type) {
    case "configuration load complete":
      const enabled =
        action.config?.parameters?.annotations_cell_ontology_enabled;
      const terms = action.config?.parameters?.annotations_cell_ontology_terms;
      return {
        ...state,
        loading: false,
        enabled,
        terms
      };
    default:
      return state;
  }
};

export default Ontology;
