import { RawSchema } from "../../../../src/common/types/entities";

export const schema: { schema: RawSchema } = {
  schema: {
    annotations: {
      obs: {
        columns: [
          {
            name: "name_0",
            type: "string",
            writable: false,
          },
          {
            name: "n_genes",
            type: "int32",
            writable: false,
          },
          {
            name: "percent_mito",
            type: "float32",
            writable: false,
          },
          {
            name: "n_counts",
            type: "float32",
            writable: false,
          },
          {
            name: "louvain",
            type: "string",
            writable: false,
          },
        ],
        index: "name_0",
      },
      var: {
        columns: [
          {
            name: "name_0",
            type: "string",
            writable: false,
          },
          {
            name: "n_cells",
            type: "int32",
            writable: false,
          },
        ],
        index: "name_0",
      },
    },
    dataframe: {
      nObs: 2638,
      nVar: 1838,
      type: "float32",
    },
    layout: {
      obs: [
        {
          dims: ["draw_graph_fr_0", "draw_graph_fr_1"],
          name: "draw_graph_fr",
          type: "float32",
        },
        {
          dims: ["pca_0", "pca_1"],
          name: "pca",
          type: "float32",
        },
        {
          dims: ["tsne_0", "tsne_1"],
          name: "tsne",
          type: "float32",
        },
        {
          dims: ["umap_0", "umap_1"],
          name: "umap",
          type: "float32",
        },
      ],
    },
  },
};
