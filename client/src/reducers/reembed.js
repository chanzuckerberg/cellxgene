/*
controller state is not part of the undo/redo history
*/
export const reembedController = (
  state = {
    pendingFetch: null,
  },
  action
) => {
  switch (action.type) {
    case "reembed: request start": {
      return {
        ...state,
        pendingFetch: action.abortableFetch,
      };
    }
    case "reembed: request aborted":
    case "reembed: request cancel":
    case "reembed: request completed": {
      return {
        ...state,
        pendingFetch: null,
      };
    }
    default: {
      return state;
    }
  }
};

const defaultPrepParams = {
  doPreprocess: false,
  minCountsCF: 0,
  minGenesCF: 0,
  minCellsGF: 0,
  maxCellsGF: 100,
  minCountsGF: 0,
  doSAM: false,
  nTopGenesHVG: 2000,
  nBinsHVG: 20,
};
const defaultBatchParams = {
  doBatch: false,
  batchMethod: "Scanorama",
  batchKey: "",
  scanoramaKnn: 20,
  scanoramaSigma: 15,
  scanoramaAlpha: 0.1,
  scanoramaBatchSize: 5000,
  bbknnNeighborsWithinBatch: 3,
};
const defaultDimredParams = {
  numPCs: 150,
  pcaSolver: "randomized",
  neighborsKnn: 20,
  neighborsMethod: "umap",
  distanceMetric: "cosine",
  nnaSAM: 50,
  weightModeSAM: "dispersion",
  umapMinDist: 0.1,
  logTransform: false,
  scaleData: false,
  dataLayer: "X",
  dataLayerExpr: "X",
  sumNormalizeCells: false,
};
export const defaultReembedParams = {
  ...defaultPrepParams,
  ...defaultBatchParams,
  ...defaultDimredParams,
  displayName: {}
};
export const reembedParameters = (state = defaultReembedParams, action) => {
  switch (action.type) {
    case "reembed: load": {
      const { params } = action;
      return params;
    }
    case "reembed: set parameter": {
      const { key, value } = action;
      return {
        ...state,
        [key]: value,
      };
    }
    case "reembed: set parameters": {
      const { params } = action;
      return {
        ...state,
        ...params,
      };
    }
    default:
      return state;
  }
};
