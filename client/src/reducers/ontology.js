const Ontology = (
  state = {
    enabled: false, // are ontology terms enabled?
    terms: null, // an array of term names, eg, ['cell', 'lung cell', ...]
    termSet: null, // a Set object containing all terms, for fast lookup
    loading: true,
  },
  action
) => {
  switch (action.type) {
    case "configuration load complete": {
      const enabled =
        action.config?.parameters?.annotations_cell_ontology_enabled ?? false;
      const terms = action.config?.parameters?.annotations_cell_ontology_terms;

      const termSet = new Set(terms);
      return {
        ...state,
        loading: false,
        enabled,
        terms,
        termSet,
      };
    }
    default: {
      return state;
    }
  }
};

export default Ontology;
