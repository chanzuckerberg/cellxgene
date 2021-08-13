export const datasets = {
  pbmc3k: {
    title: "pbmc3k",
    dataframe: {
      nObs: "2638",
      nVar: "1838",
      type: "float32",
    },
    categorical: {
      louvain: {
        "B cells": "342",
        "CD14+ Monocytes": "480",
        "CD4 T cells": "1144",
        "CD8 T cells": "316",
        "Dendritic cells": "37",
        "FCGR3A+ Monocytes": "150",
        Megakaryocytes: "15",
        "NK cells": "154",
      },
    },
    continuous: {
      n_genes: "int32",
      percent_mito: "float32",
      n_counts: "float32",
    },
    cellsets: {
      lasso: [
        {
          "coordinates-as-percent": { x1: 0.1, y1: 0.25, x2: 0.7, y2: 0.75 },
          count: "1131",
        },
      ],
      categorical: [
        {
          metadata: "louvain",
          values: ["B cells", "Megakaryocytes"],
          count: "357",
        },
      ],
      continuous: [
        {
          metadata: "n_genes",
          "coordinates-as-percent": { x1: 0.25, y1: 0.5, x2: 0.55, y2: 0.5 },
          count: "1537",
        },
      ],
    },

    diffexp: {
      cellset1: [
        { kind: "categorical", metadata: "louvain", values: ["B cells"] },
      ],
      cellset2: [
        {
          kind: "categorical",
          metadata: "louvain",
          values: ["CD4 T cells", "NK cells"],
        },
      ],
      "gene-results": [
        "HLA-DRB1",
        "HLA-DPB1",
        "CD79A",
        "HLA-DPA1",
        "HLA-DQA1",
        "CD79B",
        "HLA-DQB1",
        "MS4A1",
        "IL32",
        "CD37",
      ],
    },

    genes: {
      bulkadd: ["S100A8", "FCGR3A", "LGALS2", "GSTP1"],
      search: "ACD",
    },
    subset: {
      cellset1: [
        {
          kind: "categorical",
          metadata: "louvain",
          values: ["B cells", "Megakaryocytes"],
        },
      ],
      count: "357",
      categorical: {
        louvain: {
          "B cells": "342",
          "CD14+ Monocytes": "0",
          "CD4 T cells": "0",
          "CD8 T cells": "0",
          "Dendritic cells": "0",
          "FCGR3A+ Monocytes": "0",
          Megakaryocytes: "15",
          "NK cells": "0",
        },
      },
      lasso: {
        "coordinates-as-percent": { x1: 0.25, y1: 0.05, x2: 0.75, y2: 0.55 },
        count: "331",
      },
    },
    scatter: {
      genes: { x: "S100A8", y: "FCGR3A" },
    },
    pan: {
      "coordinates-as-percent": { x1: 0.75, y1: 0.75, x2: 0.35, y2: 0.35 },
    },
    features: {
      panzoom: {
        lasso: {
          "coordinates-as-percent": { x1: 0.3, y1: 0.3, x2: 0.5, y2: 0.5 },
          count: "38",
        },
      },
    },
    categoryLabel: {
      lasso: {
        "coordinates-as-percent": { x1: 0.05, y1: 0.3, x2: 0.5, y2: 0.5 },
      },
      newCount: {
        bySubsetConfig: {
          false: "668",
          true: "659",
        },
      },
    },
    annotationsFromFile: {
      count: {
        bySubsetConfig: {
          false: "1161",
          true: "852",
        },
      },
    },
    clip: {
      min: "30",
      max: "70",
      metadata: "n_genes",
      gene: "S100A8",
      "coordinates-as-percent": { x1: 0.25, y1: 0.5, x2: 0.55, y2: 0.5 },
      count: "386",
      "gene-cell-count": "416",
    },
  },
};
