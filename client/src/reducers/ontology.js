import TEST_ONTOLOGIES from "./TEST_ONTOLOGIES";

// jshint esversion: 6
const Ontology = (
  state = {
    data: null,
    loading: null
  },
  action
) => {
  switch (action.type) {
    case "configuration load complete":
      return {
        ...state,
        loading: false,
        data: TEST_ONTOLOGIES
      };
    default:
      return state;
  }
};

export default Ontology;
