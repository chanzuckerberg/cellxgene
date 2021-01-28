/**
 * Gene set state. Geneset UI state is in a different reducer.
 *
 * geneset reducer state is a Map object, where:
 *  key: the gene set name, a string.
 *  val: the genes in the gene set, which is a Set of string
 *
 * Ie, Map<string, Set<string>>
 */
const GeneSets = (
  state = {
    initialized: false,
    genesets: new Map(),
  },
  action
) => {
  switch (action.type) {
    /**
     * Initial, load-time bootstrap.
     * {
     *    type: "geneset: initial load"
     *    init: Array<Tuple<string, Array<string>>>
     * }
     */
    case "geneset: initial load": {
      const { init } = action;
      const genesets = new Map();
      for (const gs of init) {
        const genes = new Set();
        for (const gene of gs[1]) {
          genes.add(gene);
        }
        genesets.set(gs[0], genes);
      }
      return {
        genesets,
        initialized: true,
      };
    }

    /**
     * {
     *    type: "geneset: create",
     *    name: string, // gene set name
     *    genes: Set<string> || Array<string> || undefined
     * }
     */
    case "geneset: create": {
      const { name, genes } = action;
      if (state.genesets.has(name))
        throw new Error("geneset: create -- name already defined.");

      const genesets = new Map(state.genesets); // clone
      genesets.set(name, new Set(genes || []));
      return {
        ...state,
        genesets,
      };
    }

    /**
     * {
     *    type: "geneset: delete",
     *    name: string
     * }
     */
    case "geneset: delete": {
      const { name } = action;
      if (!state.genesets.has(name))
        throw new Error("geneset: delete -- name does not exist.");

      const genesets = new Map(state.genesets); // clone
      genesets.delete(name);
      return {
        ...state,
        genesets,
      };
    }

    /**
     * {
     *    type: "geneset: add genes"
     *    name: string, // gene set name
     *    genes: Set<string> || Array<string>  // genes to add
     * }
     */
    case "geneset: add genes": {
      const { name, genes } = action;
      if (!state.genesets.has(name))
        throw new Error("geneset: add genes -- name does not exist.");

      const genesets = new Map(state.genesets); // clone
      const newGenes = new Set(state.genesets.get(name)); // clone
      for (const gene of genes) {
        newGenes.add(gene);
      }
      genesets.set(name, newGenes);
      return {
        ...state,
        genesets,
      };
    }

    /**
     * {
     *    type: "geneset: del genes"
     *    name: string, // gene set name
     *    genes: Set<string> || Array<string>  // genes to delete
     * }
     */
    case "geneset: del genes": {
      const { name, genes } = action;

      if (!state.genesets.has(name))
        throw new Error("geneset: add genes -- name does not exist.");

      const genesets = new Map(state.genesets); // clone
      const newGenes = new Set(state.genesets.get(name)); // clone
      for (const gene of genes) {
        newGenes.delete(gene);
      }
      genesets.set(name, newGenes);
      return {
        ...state,
        genesets,
      };
    }

    default:
      return state;
  }
};

export default GeneSets;
