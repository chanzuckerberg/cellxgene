/**
 * Sample set state. Sampleset UI state is in a different reducer.
 *
 * sampleset reducer state is a Map object, where:
 *  key: the sampleset name, a string.
 *  val: the sampleset defined as an object ("geneset object")
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
 * Sampleset and samples Map order is significant, and will be preserved across
 * CRUD operations on either.
 *
 * This reducer does light error checking, but not as much as the backend
 * routes. Do not rely on it to enforce sampleset integrity - eg, no duplicate
 * samples in a sampleset.
 */
import { diffexpPopNamePrefix1, diffexpPopNamePrefix2 } from "../globals";

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
     *    type: "sampleset: initial load"
     *    data: JSON response
     * }
     */
    case "sampleset: initial load": {
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
            geneDescription: gene?.gene_description ?? "",
          });
        }
        const gs = {
          genesetName: gsData.geneset_name,
          genesetDescription: gsData?.geneset_description ?? "",
          genes,
        };
        genesets.set(gsData.geneset_name, gs);
      }

      return {
        ...state,
        initialized: true,
        lastTid,
        genesets,
      };
    }

    /**
     * Creates a new & empty sampleset with the given name and description.
     * {
     *    type: "sampleset: create",
     *    genesetName: string, // sample set name
     *    genesetDescription: string, // sampleset description
     * }
     *
     */
    case "sampleset: create": {
      const { genesetName, genesetDescription } = action;
      if (
        typeof genesetName !== "string" ||
        !genesetName ||
        genesetDescription === undefined
      )
        throw new Error(
          "sampleset: create -- name or description unspecified."
        );
      if (state.genesets.has(genesetName))
        throw new Error("sampleset: create -- name already defined.");

      const genesets = new Map([
        [
          genesetName,
          {
            genesetName,
            genesetDescription,
            genes: new Map(),
          },
        ],
        ...state.genesets,
      ]); // clone and add new sampleset to beginning

      return {
        ...state,
        genesets,
      };
    }

    /**
     * Deletes the named sampleset, if it exists. Throws if it does not.
     * {
     *    type: "sampleset: delete",
     *    genesetName: string
     * }
     */
    case "sampleset: delete": {
      const { genesetName } = action;
      if (!state.genesets.has(genesetName))
        throw new Error("sampleset: delete -- geneset name does not exist.");

      const genesets = new Map(state.genesets); // clone
      genesets.delete(genesetName);
      return {
        ...state,
        genesets,
      };
    }

    /**
     * Update the named sampleset with a new name and description.  Preserves the existing
     * order of the sampleset, even when the samplesetName changes.
     * {
     *    type: "sampleset: update",
     *    genesetName: string, current name of geneset to be updated
     *    update: {
     *        genesetName: string, new name
     *        genesetDescription: string, new description
     *    }
     * }
     *
     * For example, if you want to update JUST the description:
     *   dispatch({
     *     action: "sampleset: update",
     *     genesetName: "foo",
     *     update: { genesetName: "foo", genesetDescription: "a new description"}
     *   })
     */
    case "sampleset: update": {
      const { genesetName, update } = action;

      if (
        typeof genesetName !== "string" ||
        !genesetName ||
        !state.genesets.has(genesetName)
      )
        throw new Error(
          "sampleset: update -- geneset name unspecified or does not exist."
        );

      /* now that we've confirmed the gene set exists, check for duplicates */
      const genesetNameIsDuplicate = state.genesets.has(update.genesetName);
      const descriptionIsDuplicate =
        state.genesets.get(update.genesetName) &&
        state.genesets.get(update.genesetName).genesetDescription ===
          update.genesetDescription;

      if (genesetNameIsDuplicate && descriptionIsDuplicate)
        throw new Error(
          "sampleset: update -- update specified existing name and description."
        );

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
     * Adds samples to the sampleset.  They are appended to the END of the sampleset, in the
     * order provided. Duplicates or samples already in the sampleset, will be ignored.
     * {
     *    type: "sampleset: add samples"
     *    genesetName: <string>, // sample set name
     *    genes: Array<{
     *      geneSymbol: <string>,
     *      geneDescription: <string>
     *    }>
     * }
     *
     * Example:
     *   dispatch({
     *     type: "add genes",
     *     samplesetName: "foo",
     *     genes: [ { geneSymbol: "FOXP", geneDescription: "test" }]
     *   });
     */
    case "sampleset: add samples": {
      const { genesetName, genes } = action;

      if (!state.genesets.has(genesetName))
        throw new Error(
          "sampleset: add samples -- sampleset name does not exist."
        );

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
     * Delete samples from the named sampleset. Will throw if the samplesetName does
     * not exist.  Will ignore geneSymbols that do not exist.
     * {
     *    type: "sampleset: delete samples",
     *    genesetName: <string>, // the geneset from which to delete genes
     *    geneSymbols: [<string>, ...], // the gene symbols to delete.
     * }
     *
     * Example:
     *  dispatch({
     *    type: "sampleset: delete samples",
     *    genesetName: "a geneset name",
     *    geneSymbols: ["F5"]
     *  })
     */
    case "sampleset: delete samples": {
      const { genesetName, geneSymbols } = action;
      if (!state.genesets.has(genesetName))
        throw new Error(
          "sampleset: delete samples -- geneset name does not exist."
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
     * Set/update the description of the sample.  NOTE that this does not allow the name
     * of the sample to change - only "sampleset: add" and "sampleset: delete" can change
     * the samples in a sampleset.  Use this to update a sample description AFTER you add it
     * to the sampleset.
     * {
     *    type: "sampleset: set sample description",
     *    genesetName: <string>, // the sampleset to update
     *    update: {
     *      geneSymbol: <string>, // the sample to update, MUST exist already in the sampleset
     *      geneDescription: <string>
     *    }
     * }
     *
     * Example:
     *  dispatch({
     *      type: "sampleset: set sample description",
     *      genesetName: "my fav geneset",
     *      update: {
     *        geneSymbol: "F5",
     *        geneDescription: "tada, moar description"
     *      }
     *  })
     */
    case "sampleset: set sample description": {
      const { genesetName, update } = action;
      if (!state.genesets.has(genesetName))
        throw new Error(
          "sampleset: set sample description -- geneset name does not exist."
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
        throw new Error("sampleset: set sample description -- no such gene");
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
    case "sampleset: set tid": {
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

    case "request differential expression success": {
      const { data } = action;

      const dateString = new Date().toLocaleString();

      const genesetNames = {
        positive: `${diffexpPopNamePrefix1} (${dateString})`,
        negative: `${diffexpPopNamePrefix2} (${dateString})`,
      };

      const diffExpGeneSets = [];
      for (const polarity of Object.keys(genesetNames)) {
        const genes = new Map(
          data[polarity].map((diffExpGene) => [
            diffExpGene[0],
            {
              geneSymbol: diffExpGene[0],
            },
          ])
        );
        diffExpGeneSets.push([
          genesetNames[polarity],
          {
            genesetName: genesetNames[polarity],
            genesetDescription: "",
            genes,
          },
        ]);
      }

      const genesets = new Map([...diffExpGeneSets, ...state.genesets]); // clone

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
