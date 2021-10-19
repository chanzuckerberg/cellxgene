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

export const preprocessController = (
  state = {
    pendingFetch: null,
  },
  action
) => {
  switch (action.type) {
    case "preprocess: request start": {
      return {
        ...state,
        pendingFetch: action.abortableFetch,
      };
    }
    case "preprocess: request aborted":
    case "preprocess: request cancel":
    case "preprocess: request completed": {
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

export const defaultPrepParams = {
  doPreprocess: false,
  minCountsCF: 0,
  minGenesCF: 0,
  minCellsGF: 0,
  maxCellsGF: 100,
  minCountsGF: 0,
  nTopGenesHVG: 2000,
  nBinsHVG: 20,
  doBatchPrep: false,
  batchPrepKey: "",
  batchPrepLabel: "",
  dataLayer: "X",
  logTransform: false,
  sumNormalizeCells: false  
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
  doSAM: false,
  nnaSAM: 50,
  scaleData: false,
  weightModeSAM: "rms",
  umapMinDist: 0.1,
  dataLayerExpr: "X"
};
export const defaultReembedParams = {
  ...defaultPrepParams,
  ...defaultBatchParams,
  ...defaultDimredParams,
  batchPrepParams: {}
};
export const reembedParameters = (state = defaultReembedParams, action) => {
  switch (action.type) {
    case "reembed: load": {
      const { params } = action;
      return params;
    }
    case "reembed: set parameter": {
      const { batchPrepParams, batchPrepKey, batchPrepLabel } = state;
      const { key, value } = action;
      if (key === "doBatchPrep" && !value){
        return {
          ...state,
          batchPrepParams: {},
          [key]: value,
          batchPrepKey: "",
          batchPrepLabel: ""
        }
      }else if (key === "batchPrepKey" && value !== ""){ // create new param dict for batch key
        batchPrepParams[value] = {};
        return {
          ...state,
          [key]: value,
          batchPrepParams
        };                       
      }else if (key === "batchPrepLabel" && value !== ""){
        if (value.toString() in batchPrepParams[batchPrepKey]){
          batchPrepParams[batchPrepKey][value.toString()] = {...batchPrepParams[batchPrepKey][value.toString()],
                                                             batchPrepKey, doBatchPrep: true, batchPrepLabel: value.toString()};
        } else {
          batchPrepParams[batchPrepKey][value.toString()] = {...defaultPrepParams, batchPrepKey, doBatchPrep: true, batchPrepLabel: value.toString()};
        }
        
        return {
          ...state,
          [key]: value.toString(),
          batchPrepParams
        };                
      } else if (key !== "batchPrepLabel" && key !== "batchPrepKey" && (key in defaultPrepParams) && batchPrepLabel !== "") {
        batchPrepParams[batchPrepKey][batchPrepLabel] = {...batchPrepParams[batchPrepKey][batchPrepLabel], [key]: value}
        return {
          ...state,
          batchPrepParams
        };        
      } else if (key === "doBatch" && !value) {
        return {
          ...state,
          [key]: value,
          batchKey: ""
        }
        
      } else {
        return {
          ...state,
          [key]: value,
        };
      }

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
