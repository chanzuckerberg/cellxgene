import * as globals from "../globals";
import { AnnoMatrixLoader, AnnoMatrixObsCrossfilter } from "../annoMatrix";
import { postExplainNewTab } from "../components/framework/toasters";
import {
  KEYS,
  storageGet,
  storageSet,
  WORK_IN_PROGRESS_WARN_STATE,
} from "../components/util/localStorage";
import {
  catchErrorsWrap,
  doJsonRequest,
  dispatchNetworkErrorMessageToUser,
} from "../util/actionHelpers";
import { loadUserColorConfig } from "../util/stateManager/colorHelpers";
import * as selnActions from "./selection";
import * as annoActions from "./annotation";
import * as viewActions from "./viewStack";
import * as embActions from "./embedding";
import {
  createDatasetUrl,
  createExplorerUrl,
  createAPIPrefix,
} from "../util/stateManager/collectionsHelpers";
import * as genesetActions from "./geneset";

// eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
function setGlobalConfig(config: any) {
  /**
   * Set any global run-time config not _exclusively_ managed by the config reducer.
   * This should only set fields defined in globals.globalConfig.
   */
  globals.globalConfig.maxCategoricalOptionsToDisplay =
    config?.parameters?.["max-category-items"] ??
    globals.globalConfig.maxCategoricalOptionsToDisplay;
}

/*
return promise fetching user-configured colors
*/
// eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
async function userColorsFetchAndLoad(dispatch: any) {
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

// eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
async function configFetch(dispatch: any) {
  return fetchJson("config").then((response) => {
    const config = { ...globals.configDefaults, ...response.config };

    setGlobalConfig(config);

    dispatch({
      type: "configuration load complete",
      config,
    });
    return config;
  });
}

// eslint-disable-next-line @typescript-eslint/no-explicit-any -- - FIXME: disabled temporarily on migrate to TS.
async function collectionFetchAndLoad(dispatch: any) {
  /*
  Fetch dataset meta for the current visualization then fetch the corresponding collection.
   */
  const datasetMeta = await datasetMetaFetch();
  const { collection_id: collectionId, dataset_id: selectedDatasetId } =
    datasetMeta;
  const collection = await collectionFetch(collectionId);
  dispatch({
    type: "collection load complete",
    collection,
    selectedDatasetId,
  });
}

// eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
async function userInfoFetch(dispatch: any) {
  return fetchJson("userinfo").then((response) => {
    const { userinfo: userInfo } = response || {};
    dispatch({
      type: "userInfo load complete",
      userInfo,
    });
    return userInfo;
  });
}

// eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
async function genesetsFetch(dispatch: any, config: any) {
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

async function datasetMetaFetch() {
  /*
   Fetch dataset meta for the current dataset.
   TODO(cc) revisit swap of explorer URL origin for environments without a corresponding Portal instance (eg local, canary)
   */
  const explorerUrl = createExplorerUrl();
  const explorerUrlParam = encodeURIComponent(explorerUrl);
  return fetchPortalJson(`datasets/meta?url=${explorerUrlParam}`);
}

async function collectionFetch(collectionId: string) {
  /*
   Fetch collection with the given ID.
   */
  return fetchPortalJson(`collections/${collectionId}`);
}

// eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
function prefetchEmbeddings(annoMatrix: any) {
  /*
  prefetch requests for all embeddings
  */
  const { schema } = annoMatrix;
  // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
  const available = schema.layout.obs.map((v: any) => v.name);
  // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
  available.forEach((embName: any) => annoMatrix.prefetch("emb", embName));
}

/*
Application bootstrap
*/
// eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types --- FIXME: disabled temporarily on migrate to TS.
const doInitialDataLoad = () =>
  // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
  catchErrorsWrap(async (dispatch: any) => {
    dispatch({ type: "initial data load start" });

    try {
      const [config, schema] = await Promise.all([
        configFetch(dispatch),
        // @ts-expect-error ts-migrate(2554) FIXME: Expected 0 arguments, but got 1.
        schemaFetch(dispatch),
        userColorsFetchAndLoad(dispatch),
        userInfoFetch(dispatch),
        collectionFetchAndLoad(dispatch),
      ]);

      genesetsFetch(dispatch, config);

      const baseDataUrl = `${globals.API.prefix}${globals.API.version}`;
      const annoMatrix = new AnnoMatrixLoader(baseDataUrl, schema.schema);
      const obsCrossfilter = new AnnoMatrixObsCrossfilter(annoMatrix);
      prefetchEmbeddings(annoMatrix);

      dispatch({
        type: "annoMatrix: init complete",
        annoMatrix,
        obsCrossfilter,
      });
      dispatch({ type: "initial data load complete" });

      const defaultEmbedding = config?.parameters?.default_embedding;
      const layoutSchema = schema?.schema?.layout?.obs ?? [];
      if (
        defaultEmbedding &&
        // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
        layoutSchema.some((s: any) => s.name === defaultEmbedding)
      ) {
        dispatch(embActions.layoutChoiceAction(defaultEmbedding));
      }
    } catch (error) {
      dispatch({ type: "initial data load error", error });
    }
  }, true);

// eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any -- - FIXME: disabled temporarily on migrate to TS.
function requestSingleGeneExpressionCountsForColoringPOST(gene: any) {
  return {
    type: "color by expression",
    gene,
  };
}

// eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any -- - FIXME: disabled temporarily on migrate to TS.
const requestUserDefinedGene = (gene: any) => ({
  type: "request user defined gene success",

  data: {
    genes: [gene],
  },
});

// eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
const dispatchDiffExpErrors = (dispatch: any, response: any) => {
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

const requestDifferentialExpression =
  (
    // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any -- - FIXME: disabled temporarily on migrate to TS.
    set1: any,
    // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any -- - FIXME: disabled temporarily on migrate to TS.
    set2: any,
    num_genes = 50
    // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any -- - FIXME: disabled temporarily on migrate to TS.
  ) =>
  async (dispatch: any, getState: any) => {
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
        // @ts-expect-error ts-migrate(7053) FIXME: Element implicitly has an 'any' type because expre... Remove this comment to see the full error message
        // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
        diffexpLists[polarity] = response[polarity].map((v: any) => [
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

// eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any -- - FIXME: disabled temporarily on migrate to TS.
export const checkExplainNewTab = () => (dispatch: any) => {
  /*
  Opens toast "work in progress" warning.
   */
  if (
    storageGet(KEYS.WORK_IN_PROGRESS_WARN) === WORK_IN_PROGRESS_WARN_STATE.ON
  ) {
    dispatch({ type: "work in progress warning displayed" });
    postExplainNewTab(
      "To maintain your in-progress work on the previous dataset, we opened this dataset in a new tab."
    );
    storageSet(KEYS.WORK_IN_PROGRESS_WARN, WORK_IN_PROGRESS_WARN_STATE.OFF);
  }
};

// eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any -- - FIXME: disabled temporarily on migrate to TS.
export const openDataset = (dataset: any) => (dispatch: any) => {
  /*
  Update in a new tab the browser location to dataset's deployment URL, kick off data load.
  */

  const deploymentUrl = dataset.dataset_deployments?.[0].url ?? "";
  const datasetUrl = createDatasetUrl(deploymentUrl);

  dispatch({ type: "dataset opened" });
  storageSet(KEYS.WORK_IN_PROGRESS_WARN, WORK_IN_PROGRESS_WARN_STATE.ON);
  window.open(datasetUrl, "_blank");
};

// eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any -- - FIXME: disabled temporarily on migrate to TS.
export const switchDataset = (dataset: any) => (dispatch: any) => {
  /*
  Update browser location to dataset's deployment URL, kick off data load. 
  TODO(cc) revisit:
    - origin (and data root) switch for environments without corresponding Portal instance (eg local, canary)
    - globals update: move to server-side, split from initial doc returned from server?
   */
  dispatch({ type: "dataset switch" });

  const deploymentUrl = dataset.dataset_deployments?.[0].url ?? "";
  const datasetUrl = createDatasetUrl(deploymentUrl);
  dispatch(updateLocation(datasetUrl));

  globals.API.prefix = createAPIPrefix(globals.API.prefix, datasetUrl);
  dispatch(doInitialDataLoad());
};

// eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any -- - FIXME: disabled temporarily on migrate to TS.
const updateLocation = (url: string) => (dispatch: any) => {
  /*
  Add entry to the session's history stack.
   */
  dispatch({ type: "location update" });
  window.history.pushState(null, "", url);
};

// eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
function fetchJson(pathAndQuery: any) {
  return doJsonRequest(
    `${globals.API.prefix}${globals.API.version}${pathAndQuery}`
  );
}

function fetchPortalJson(url: string) {
  /* 
  Fetch JSON from Portal API.
  TODO(cc) revisit - required for dataset meta and collection requests from Portal 
   */
  return doJsonRequest(`${globals.API.portalPrefix}${url}`);
}

export default {
  doInitialDataLoad,
  requestDifferentialExpression,
  requestSingleGeneExpressionCountsForColoringPOST,
  requestUserDefinedGene,
  checkExplainNewTab,
  openDataset,
  switchDataset,
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
  annotationDeleteLabelFromCategory:
    annoActions.annotationDeleteLabelFromCategory,
  annotationRenameLabelInCategory: annoActions.annotationRenameLabelInCategory,
  annotationLabelCurrentSelection: annoActions.annotationLabelCurrentSelection,
  saveObsAnnotationsAction: annoActions.saveObsAnnotationsAction,
  saveGenesetsAction: annoActions.saveGenesetsAction,
  needToSaveObsAnnotations: annoActions.needToSaveObsAnnotations,
  layoutChoiceAction: embActions.layoutChoiceAction,
  setCellSetFromSelection: selnActions.setCellSetFromSelection,
  genesetDelete: genesetActions.genesetDelete,
  genesetAddGenes: genesetActions.genesetAddGenes,
  genesetDeleteGenes: genesetActions.genesetDeleteGenes,
};
