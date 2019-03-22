export const datasets = {
  pbmc3k: {
    title: "cellxgene: pbmc3k",
    dataframe: {
      nObs: "2638",
      nVar: "1838",
      type: "float32"
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
        "NK cells": "154"
      }
    },
    continuous: {
      n_genes: "int32",
      percent_mito: "float32",
      n_counts: "float32"
    },
    cellsets: {
      lasso: [
        {
          "coordinates-as-percent": { x1: 0.25, y1: 0.25, x2: 0.35, y2: 0.35 },
          count: "26"
        }
      ],
      categorical: [
        {
          metadata: "louvain",
          values: ["B cells", "Megakaryocytes"],
          count: "357"
        }
      ],
      continuous: [
        {
          metadata: "n_genes",
          "coordinates-as-percent": { x1: 0.25, y1: 0.5, x2: 0.55, y2: 0.5 },
          count: "1537"
        }
      ]
    },

    diffexp: {
      cellset1: [
        { kind: "categorical", metadata: "louvain", values: ["B cells"] }
      ],
      cellset2: [
        {
          kind: "categorical",
          metadata: "louvain",
          values: ["CD4 T cells", "NK cells"]
        }
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
        "CD37"
      ]
    },

    genes: {
      "bulk add": ["S100A8", "FCGR3A", "LGALS2", "GSTP1"],
      search: "ACD"
    }
  }
};
