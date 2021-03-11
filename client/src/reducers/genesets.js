/**
 * Gene set state. Geneset UI state is in a different reducer.
 *
 * geneset reducer state is a Map object, where:
 *  key: the geneset name, a string.
 *  val: the geneset defined as an object ("geneset object")
 *
 * A geneset object is:
 * {
 *    genesetName: <string>  # same as the map key
 *    genesetDescription: <string>
 *    genes: Map<<string>, {
 *      geneSymbol: <string>,  # same as the map key
 *      geneDescription: <string>
 *    }>
 * }
 *
 * Geneset and genes Map order is significant, and will be preserved across
 * CRUD operations on either.
 *
 * This reducer does light error checking, but not as much as the backend
 * routes. Do not rely on it to enforce geneset integrity - eg, no duplicate
 * genes in a geneset.
 */
const GeneSets = (
  state = {
    initialized: false,
    lastTid: undefined,
    genesets: new Map(),
  },
  action
) => {
  switch (action.type) {
    /**
     * Initial, load-time bootstrap.
     * {
     *    type: "geneset: initial load"
     *    data: JSON response
     * }
     */
    case "geneset: initial load": {
      const { data } = action;

      if (
        !data ||
        typeof data.tid !== "number" ||
        !Array.isArray(data.genesets)
      )
        throw new Error("missing or malformed JSON response");

      const lastTid = data.tid;
      const genesetsData = data.genesets;
      const genesets = new Map();

      for (const gsData of genesetsData) {
        const genes = new Map();
        for (const gene of gsData.genes) {
          genes.set(gene.gene_symbol, {
            geneSymbol: gene.gene_symbol,
            geneDescription: gene?.["gene_description"] ?? "",
          });
        }
        const gs = {
          genesetName: gsData.geneset_name,
          genesetDescription: gsData?.["geneset_description"] ?? "",
          genes,
        };
        genesets.set(gsData.geneset_name, gs);
      }

      return {
        initialized: true,
        lastTid,
        genesets,
      };
    }

    /**
     * Creates a new & empty geneset with the given name and description.
     * {
     *    type: "geneset: create",
     *    genesetName: string, // gene set name
     *    genesetDescription: string, // geneset description
     * }
     *
     */
    case "geneset: create": {
      const {
        genesetName,
        genesetDescription = "No description provided",
      } = action;
      if (
        typeof genesetName !== "string" ||
        !genesetName ||
        genesetDescription === undefined
      )
        throw new Error("geneset: create -- name or description unspecified.");
      if (state.genesets.has(genesetName))
        throw new Error("geneset: create -- name already defined.");

      const genesets = new Map(state.genesets); // clone
      genesets.set(genesetName, {
        genesetName,
        genesetDescription,
        genes: new Map(),
      });

      return {
        ...state,
        genesets,
      };
    }

    /**
     * Deletes the named geneset, if it exists. Throws if it does not.
     * {
     *    type: "geneset: delete",
     *    genesetName: string
     * }
     */
    case "geneset: delete": {
      const { genesetName } = action;
      if (!state.genesets.has(genesetName))
        throw new Error("geneset: delete -- geneset name does not exist.");

      const genesets = new Map(state.genesets); // clone
      genesets.delete(genesetName);
      return {
        ...state,
        genesets,
      };
    }

    /**
     * Update the named geneset with a new name and description.  Preserves the existing
     * order of the geneset, even when the genesetName changes.
     * {
     *    type: "geneset: update",
     *    genesetName: string, current name of geneset to be updated
     *    update: {
     *        genesetName: string, new name
     *        genesetDescription: string, new description
     *    }
     * }
     *
     * For example, if you want to update JUST the description:
     *   dispatch({
     *     action: "geneset: update",
     *     genesetName: "foo",
     *     update: { genesetName: "foo", genesetDescription: "a new description"}
     *   })
     */
    case "geneset: update": {
      const { genesetName, update } = action;
      if (
        typeof genesetName !== "string" ||
        !genesetName ||
        !state.genesets.has(genesetName)
      )
        throw new Error(
          "geneset: update -- geneset name unspecified or does not exist."
        );
      if (state.genesets.has(update.genesetName))
        throw new Error("geneset: update -- update specified existing name.");

      const prevGs = state.genesets.get(genesetName);
      const newGs = {
        ...update,
        genes: prevGs.genes,
      }; // clone

      // clone the map, preserving current insert order, but mapping name->newName.
      const genesets = new Map();
      for (const [name, gs] of state.genesets) {
        if (name === genesetName) genesets.set(newGs.genesetName, newGs);
        else genesets.set(name, gs);
      }

      return {
        ...state,
        genesets,
      };
    }

    /**
     * Adds genes to the geneset.  They are appended to the END of the geneset, in the
     * order provided. Duplicates or genes already in the geneset, will be ignored.
     * {
     *    type: "geneset: add genes"
     *    genesetName: <string>, // gene set name
     *    genes: Array<{
     *      geneSymbol: <string>,
     *      geneDescription: <string>
     *    }>
     * }
     *
     * Example:
     *   dispatch({
     *     type: "add genes",
     *     genesetName: "foo",
     *     genes: [ { geneSymbol: "FOXP", geneDescription: "test" }]
     *   });
     */
    case "geneset: add genes": {
      const { genesetName, genes } = action;

      if (!state.genesets.has(genesetName))
        throw new Error("geneset: add genes -- geneset name does not exist.");

      // clone
      const genesets = new Map(state.genesets);
      const gs = {
        ...genesets.get(genesetName),
        genes: new Map(genesets.get(genesetName).genes),
      };
      genesets.set(genesetName, gs);

      // add
      const newGenes = gs.genes;
      for (const gene of genes) {
        const { geneSymbol } = gene;
        const geneDescription = gene?.geneDescription ?? "";
        // ignore genes already present
        if (!newGenes.has(geneSymbol))
          newGenes.set(geneSymbol, {
            geneSymbol,
            geneDescription,
          });
      }

      return {
        ...state,
        genesets,
      };
    }

    /**
     * Delete genes from the named geneset. Will throw if the genesetName does
     * not exist.  Will ignore geneSymbols that do not exist.
     * {
     *    type: "geneset: delete genes",
     *    genesetName: <string>, // the geneset from which to delete genes
     *    geneSymbols: [<string>, ...], // the gene symbols to delete.
     * }
     *
     * Example:
     *  dispatch({
     *    type: "geneset: delete genes",
     *    genesetName: "a geneset name",
     *    geneSymbols: ["F5"]
     *  })
     */
    case "geneset: delete genes": {
      const { genesetName, geneSymbols } = action;
      if (!state.genesets.has(genesetName))
        throw new Error(
          "geneset: delete genes -- geneset name does not exist."
        );

      // clone
      const genesets = new Map(state.genesets);
      const gs = {
        ...genesets.get(genesetName),
        genes: new Map(genesets.get(genesetName).genes),
      };
      genesets.set(genesetName, gs);

      // delete
      const { genes } = gs;
      for (const geneSymbol of geneSymbols) {
        genes.delete(geneSymbol);
      }
      return {
        ...state,
        genesets,
      };
    }

    /**
     * Set/update the description of the gene.  NOTE that this does not allow the name
     * of the gene to change - only "geneset: add" and "geneset: delete" can change
     * the genes in a geneset.  Use this to update a gene description AFTER you add it
     * to the geneset.
     * {
     *    type: "geneset: set gene description",
     *    genesetName: <string>, // the geneset to update
     *    update: {
     *      geneSymbol: <string>, // the gene to update, MUST exist already in the geneset
     *      geneDescription: <string>
     *    }
     * }
     *
     * Example:
     *  dispatch({
     *      type: "geneset: set gene description",
     *      genesetName: "my fav geneset",
     *      update: {
     *        geneSymbol: "F5",
     *        geneDescription: "tada, moar description"
     *      }
     *  })
     */
    case "geneset: set gene description": {
      const { genesetName, update } = action;
      if (!state.genesets.has(genesetName))
        throw new Error(
          "geneset: set gene description -- geneset name does not exist."
        );

      // clone
      const genesets = new Map(state.genesets);
      const gs = {
        ...genesets.get(genesetName),
        genes: new Map(genesets.get(genesetName).genes),
      };
      genesets.set(genesetName, gs);

      const { geneSymbol, geneDescription } = update;
      const gene = gs.genes.get(geneSymbol);
      if (!gene)
        throw new Error("geneset: set gene description -- no such gene");
      gs.genes.set(geneSymbol, {
        geneSymbol,
        geneDescription,
      });

      return {
        ...state,
        genesets,
      };
    }

    /**
     * Used by autosave to update the server synchronization TID
     */
    case "geneset: set tid": {
      const { tid } = action;
      if (!Number.isInteger(tid) || tid < 0)
        throw new Error("TID must be a positive integer number");
      if (state.lastTid !== undefined && tid < state.lastTid)
        throw new Error("TID may not be decremented.");
      return {
        ...state,
        lastTid: tid,
      };
    }

    default:
      return state;
  }
};

export default GeneSets;
