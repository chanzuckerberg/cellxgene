export const datasets = {
  pbmc3k: {
    title: "cellxgene: pbmc3k",
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
    lasso: {
      "coordinates-as-percent": { x1: 0.25, y1: 0.25, x2: 0.35, y2: 0.35 },
      "cell count": 26
    },
    diffexp: {
      cellset1: [{ category: "louvain", values: ["B cells"] }],
      cellset2: [{ category: "louvain", values: ["CD4 T cells", "NK cells"] }],
      "gene-results": [
        "HLA-DPA1",
        "HLA-DQA1",
        "HLA-DRB1",
        "HLA-DMA",
        "CST3",
        "HLA-DPB1",
        "HLA-DQB1",
        "LGALS2",
        "FCER1A",
        "LTB"
      ]
    },
    "historgram-brush": [
      {
        "coordinates-as-percent": { x1: 0.25, y1: 0.5, x2: 0.55, y2: 0.5 },
        "cell count": 1537
      }
    ],
    genes: {
      "bulk add": ["S100A8", "FCGR3A", "LGALS2", "GSTP1"],
      search: "ACD"
    }
  }
};
