import * as globals from "../globals";
import { API } from "../globals";
import { AnnoMatrixLoader, AnnoMatrixObsCrossfilter } from "../annoMatrix";
import {
  catchErrorsWrap,
  doJsonRequest,
  dispatchNetworkErrorMessageToUser,
} from "../util/actionHelpers";
import {
  requestReembed, requestPreprocessing
} from "./reembed";
import {
  requestSankey
} from "./sankey";
import {
  requestLeiden
} from "./leiden";
import {
  postNetworkErrorToast,
  postAsyncSuccessToast,
  postAsyncFailureToast,
} from "../components/framework/toasters";
import { loadUserColorConfig } from "../util/stateManager/colorHelpers";
import * as selnActions from "./selection";
import * as annoActions from "./annotation";
import * as viewActions from "./viewStack";
import * as embActions from "./embedding";
import * as genesetActions from "./geneset";
import { defaultReembedParams } from "../reducers/reembed";
import { _switchEmbedding } from "./embedding";
/*
return promise fetching user-configured colors
*/
async function userColorsFetchAndLoad(dispatch) {
  return fetchJson("colors").then((response) =>
    dispatch({
      type: "universe: user color load success",
      userColors: loadUserColorConfig(response),
    })
  );
}

async function schemaFetch() {
  return fetchJson("schema");
}

async function configFetch(dispatch) {
  return fetchJson("config").then((response) => {
    const config = { ...globals.configDefaults, ...response.config };
    dispatch({
      type: "configuration load complete",
      config,
    });
    return config;
  });
}

export async function userInfoFetch(dispatch) {
  return fetchJson("userinfo").then((response) => {
    const { userinfo: userInfo } = response || {};
    dispatch({
      type: "userInfo load complete",
      userInfo,
    });
    return userInfo;
  });
}

async function genesetsFetch(dispatch, config) {
  /* request genesets ONLY if the backend supports the feature */
  const defaultResponse = {
    genesets: [],
    tid: 0,
  };
  if (config?.parameters?.annotations_genesets ?? false) {
    fetchJson("genesets").then((response) => {
      dispatch({
        type: "geneset: initial load",
        data: response ?? defaultResponse,
      });
    });
  } else {
    dispatch({
      type: "geneset: initial load",
      data: defaultResponse,
    });
  }
}
export async function reembedParamsFetch(dispatch) {
  /* request reembedding parameters ONLY if the backend supports the feature */
  const defaultResponse = {
    reembedParams: defaultReembedParams,
  };
  try {
    fetchJson("reembed-parameters").then((response) => {
      const isEmpty = Object.keys(response.reembedParams).length === 0;
      dispatch({
        type: "reembed: load",
        params: isEmpty
          ? defaultResponse.reembedParams
          : response.reembedParams ?? defaultResponse.reembedParams,
      });
    });
  } catch (e) {
    dispatch({
      type: "reembed: load",
      data: defaultResponse,
    });
  }
}

function prefetchEmbeddings(annoMatrix) {
  /*
  prefetch requests for all embeddings
  */
  const { schema } = annoMatrix;
  const available = schema.layout.obs.map((v) => v.name);
  available.forEach((embName) => annoMatrix.prefetch("emb", embName));
}

function abortableFetch(request, opts, timeout = 0) {
  const controller = new AbortController();
  const { signal } = controller;

  return {
    abort: () => controller.abort(),
    isAborted: () => signal.aborted,
    ready: () => {
      if (timeout) {
        setTimeout(() => controller.abort(), timeout);
      }
      return fetch(request, { ...opts, signal });
    },
  };
}
export const requestSaveAnndataToFile = (saveName) => async (
  dispatch,
  getState
) => {
  try{
    const state = getState();
    const { annoMatrix, layoutChoice } = state;
    
    let cells = annoMatrix.rowIndex.labels();  
    cells = Array.isArray(cells) ? cells : Array.from(cells);

    const annos = []
    const annoNames = []
    
    for (const item of annoMatrix.schema.annotations?.obs?.columns) {
      if(item?.categories){
        let labels = await annoMatrix.fetch("obs",item.name)
        annos.push(labels)
        annoNames.push(item.name)
      }
    }
    const af = abortableFetch(
      `${API.prefix}${API.version}output`,
      {
        method: "PUT",
        headers: new Headers({
          Accept: "application/octet-stream",
          "Content-Type": "application/json",
        }),
        body: JSON.stringify({
          saveName: saveName,
          labelNames: annoNames,
          labels: annos,
          currentLayout: layoutChoice.current,
          filter: { obs: { index: cells } }
        }),
        credentials: "include",
      },
      60000
    );
    dispatch({
      type: "output data: request start",
      abortableFetch: af,
    });
    const res = await af.ready();
    postAsyncSuccessToast("Data has been successfully saved.");
    dispatch({
      type: "output data: request completed",
    });
    if (res.ok && res.headers.get("Content-Type").includes("application/json")) {      
      return true;
    }

    // else an error
    let msg = `Unexpected HTTP response ${res.status}, ${res.statusText}`;
    const body = await res.text();
    if (body && body.length > 0) {
      msg = `${msg} -- ${body}`;
    }
    throw new Error(msg);
  } catch (error) {
    dispatch({
      type: "ouput data: request aborted",
    });
    if (error.name === "AbortError") {
      postAsyncFailureToast("Data output was aborted.");
    } else {
      postNetworkErrorToast(`Data output: ${error.message}`);
    }
  }
}
export function requestDataLayerChange(dataLayer) {
  return async (_dispatch, _getState) => {
    const res = await fetch(
      `${API.prefix}${API.version}layer`,
      {
        method: "PUT",
        headers: new Headers({
          Accept: "application/octet-stream",
          "Content-Type": "application/json",
        }),
        body: JSON.stringify({
          dataLayer: dataLayer,
        }),
        credentials: "include",
        }
    );
    if (res.ok && res.headers.get("Content-Type").includes("application/json")) {
      return res;
    }
  }
}
export function requestReloadBackend() {
  return async (dispatch, getState) => {
    try{
      const { layoutChoice } = getState()
      const af = abortableFetch(
        `${API.prefix}${API.version}reload`,
        {
          method: "PUT",
          headers: new Headers({
            Accept: "application/octet-stream",
            "Content-Type": "application/json",
          }),
          body: JSON.stringify({
            currentLayout: layoutChoice.current,
          }),        
          credentials: "include",
          },
          60000
      );
      dispatch({
        type: "output data: request start",
        abortableFetch: af,
      });
      const res = await af.ready();      
      
      dispatch({
        type: "app: refresh"
      })

      dispatch({
        type: "output data: request completed",
      });

      postAsyncSuccessToast("Data has successfuly overwritten the backend.");

      if (res.ok && res.headers.get("Content-Type").includes("application/json")) {      
        return true;
      }

      // else an error
      let msg = `Unexpected HTTP response ${res.status}, ${res.statusText}`;
      const body = await res.text();
      if (body && body.length > 0) {
        msg = `${msg} -- ${body}`;
      }
      throw new Error(msg);
    } catch (error) {
      dispatch({
        type: "ouput data: request aborted",
      });
      if (error.name === "AbortError") {
        postAsyncFailureToast("Data output was aborted.");
      } else {
        postNetworkErrorToast(`Data output: ${error.message}`);
      }
    }
  }
}
/*
Application bootstrap
*/
const doInitialDataLoad = () =>
  catchErrorsWrap(async (dispatch) => {
    dispatch({ type: "initial data load start" });

    try {
      const [config, schema] = await Promise.all([
        configFetch(dispatch),
        schemaFetch(dispatch),
        userColorsFetchAndLoad(dispatch),
        userInfoFetch(dispatch),
      ]);
      genesetsFetch(dispatch, config);
      reembedParamsFetch(dispatch);
      await dispatch(requestDataLayerChange("X"));
      const baseDataUrl = `${globals.API.prefix}${globals.API.version}`;
      const annoMatrix = new AnnoMatrixLoader(baseDataUrl, schema.schema);
      const obsCrossfilter = new AnnoMatrixObsCrossfilter(annoMatrix);
      prefetchEmbeddings(annoMatrix);
      
      const layoutSchema = schema?.schema?.layout?.obs ?? [];
      if(layoutSchema.length > 0){
        const name = layoutSchema[0].name
        const base = annoMatrix.base();

        const [annoMatrixNew, obsCrossfilterNew] = await _switchEmbedding(
          base,
          obsCrossfilter,
          name,
          name
        ); 
        dispatch({
          type: "annoMatrix: init complete",
          annoMatrix: annoMatrixNew,
          obsCrossfilter: obsCrossfilterNew,
        });        
        
      } else { 
        dispatch({
          type: "annoMatrix: init complete",
          annoMatrix,
          obsCrossfilter,
        });
      }

      dispatch({ type: "initial data load complete" });

      const defaultEmbedding = config?.parameters?.default_embedding;
      if (
        defaultEmbedding &&
        layoutSchema.some((s) => s.name === defaultEmbedding)
      ) {
        dispatch(embActions.layoutChoiceAction(defaultEmbedding));
      }

 
    } catch (error) {
      dispatch({ type: "initial data load error", error });
    }
  }, true);

function requestSingleGeneExpressionCountsForColoringPOST(gene) {
  return {
    type: "color by expression",
    gene,
  };
}

const requestUserDefinedGene = (gene) => ({
  type: "request user defined gene success",
  data: {
    genes: [gene],
  },
});

const dispatchDiffExpErrors = (dispatch, response) => {
  switch (response.status) {
    case 403:
      dispatchNetworkErrorMessageToUser(
        "Too many cells selected for differential experesion calculation - please make a smaller selection."
      );
      break;
    case 501:
      dispatchNetworkErrorMessageToUser(
        "Differential expression is not implemented."
      );
      break;
    default: {
      const msg = `Unexpected differential expression HTTP response ${response.status}, ${response.statusText}`;
      dispatchNetworkErrorMessageToUser(msg);
      dispatch({
        type: "request differential expression error",
        error: new Error(msg),
      });
    }
  }
};

const requestDifferentialExpression = (set1, set2, num_genes = 50) => async (
  dispatch,
  getState
) => {
  dispatch({ type: "request differential expression started" });
  try {
    /*
    Steps:
    1. get the most differentially expressed genes
    2. get expression data for each
    */
    const { annoMatrix } = getState();
    const varIndexName = annoMatrix.schema.annotations.var.index;

    // Legal values are null, Array or TypedArray.  Null is initial state.
    if (!set1) set1 = [];
    if (!set2) set2 = [];

    // These lines ensure that we convert any TypedArray to an Array.
    // This is necessary because JSON.stringify() does some very strange
    // things with TypedArrays (they are marshalled to JSON objects, rather
    // than being marshalled as a JSON array).
    set1 = Array.isArray(set1) ? set1 : Array.from(set1);
    set2 = Array.isArray(set2) ? set2 : Array.from(set2);

    const res = await fetch(
      `${globals.API.prefix}${globals.API.version}diffexp/obs`,
      {
        method: "POST",
        headers: new Headers({
          Accept: "application/json",
          "Content-Type": "application/json",
        }),
        body: JSON.stringify({
          mode: "topN",
          count: num_genes,
          set1: { filter: { obs: { index: set1 } } },
          set2: { filter: { obs: { index: set2 } } },
        }),
        credentials: "include",
      }
    );

    if (!res.ok || res.headers.get("Content-Type") !== "application/json") {
      return dispatchDiffExpErrors(dispatch, res);
    }

    const response = await res.json();
    const varIndex = await annoMatrix.fetch("var", varIndexName);
    const diffexpLists = { negative: [], positive: [] };
    for (const polarity of Object.keys(diffexpLists)) {
      diffexpLists[polarity] = response[polarity].map((v) => [
        varIndex.at(v[0], varIndexName),
        ...v.slice(1),
      ]);
    }

    /* then send the success case action through */
    return dispatch({
      type: "request differential expression success",
      data: diffexpLists,
    });
  } catch (error) {
    return dispatch({
      type: "request differential expression error",
      error,
    });
  }
};

const selectAll = () => async (dispatch, getState) => {
  dispatch({ type: "select all observations" });
  try {
    const { obsCrossfilter: prevObsCrossfilter } = getState();
    const obsCrossfilter = await prevObsCrossfilter.selectAll();
    return dispatch({
      type: "selected all observations",
      obsCrossfilter,
    });
  } catch (error) {
    return dispatch({
      type: "error selecting all observations",
      error,
    });
  }
};

function fetchJson(pathAndQuery) {
  return doJsonRequest(
    `${globals.API.prefix}${globals.API.version}${pathAndQuery}`
  );
}

export default {
  doInitialDataLoad,
  requestDataLayerChange,
  requestReloadBackend,
  selectAll,
  requestDifferentialExpression,
  requestSingleGeneExpressionCountsForColoringPOST,
  requestUserDefinedGene,
  requestReembed,
  requestPreprocessing,
  requestSankey,
  requestLeiden,
  requestSaveAnndataToFile,
  setCellsFromSelectionAndInverseAction:
    selnActions.setCellsFromSelectionAndInverseAction,
  selectContinuousMetadataAction: selnActions.selectContinuousMetadataAction,
  selectCategoricalMetadataAction: selnActions.selectCategoricalMetadataAction,
  selectCategoricalAllMetadataAction:
    selnActions.selectCategoricalAllMetadataAction,
  graphBrushStartAction: selnActions.graphBrushStartAction,
  graphBrushChangeAction: selnActions.graphBrushChangeAction,
  graphBrushDeselectAction: selnActions.graphBrushDeselectAction,
  graphBrushCancelAction: selnActions.graphBrushCancelAction,
  graphBrushEndAction: selnActions.graphBrushEndAction,
  graphLassoStartAction: selnActions.graphLassoStartAction,
  graphLassoEndAction: selnActions.graphLassoEndAction,
  graphLassoCancelAction: selnActions.graphLassoCancelAction,
  graphLassoDeselectAction: selnActions.graphLassoDeselectAction,
  clipAction: viewActions.clipAction,
  subsetAction: viewActions.subsetAction,
  resetSubsetAction: viewActions.resetSubsetAction,
  annotationCreateCategoryAction: annoActions.annotationCreateCategoryAction,
  annotationRenameCategoryAction: annoActions.annotationRenameCategoryAction,
  annotationDeleteCategoryAction: annoActions.annotationDeleteCategoryAction,
  annotationCreateLabelInCategory: annoActions.annotationCreateLabelInCategory,
  requestFuseLabels: annoActions.requestFuseLabels,
  requestDeleteLabels: annoActions.requestDeleteLabels,
  annotationDeleteLabelFromCategory:
    annoActions.annotationDeleteLabelFromCategory,
  annotationRenameLabelInCategory: annoActions.annotationRenameLabelInCategory,
  annotationLabelCurrentSelection: annoActions.annotationLabelCurrentSelection,
  saveObsAnnotationsAction: annoActions.saveObsAnnotationsAction,
  saveGenesetsAction: annoActions.saveGenesetsAction,
  saveReembedParametersAction: annoActions.saveReembedParametersAction,
  needToSaveObsAnnotations: annoActions.needToSaveObsAnnotations,
  layoutChoiceAction: embActions.layoutChoiceAction,
  requestDeleteEmbedding: embActions.requestDeleteEmbedding,
  requestRenameEmbedding: embActions.requestRenameEmbedding,
  setCellSetFromSelection: selnActions.setCellSetFromSelection,
  setCellSetFromInputArray: selnActions.setCellSetFromInputArray,
  genesetDelete: genesetActions.genesetDelete,
  genesetAddGenes: genesetActions.genesetAddGenes,
  genesetDeleteGenes: genesetActions.genesetDeleteGenes,
};
