/*
Helper functions for querying, binding and sorting Portal dataset meta and collections.
*/

/* app dependencies */
import * as globals from "../../globals";

export function createAPIPrefix(existingPrefix, replaceWithPrefix) {
  /*
  When selecting a dataset from the dataset selector, the globals API prefix value must be updated to match the selected
  dataset's deployment URL. For example, update protocol://origin/dataRoot/current-dataset.cxg/api/v2 to
  protocol://origin/dataRoot/selected-dataset.cxg/api/v2.
  TODO(cc) revisit updating globals API prefix locally, possibly move to back end?
   */
  const newDataRootAndDeploymentId = bindDataRootAndDeploymentId(
    replaceWithPrefix
  );
  return existingPrefix.replace(
    /([a-z0-9_.-]+\/[a-z0-9_.-]+\.cxg)/i,
    newDataRootAndDeploymentId
  );
}

export function createDatasetUrl(deploymentUrl) {
  /*
  For local environments, switch out the data root of the given URL to /d/. For example, update 
  protocol://origin/dataRoot/dataset.cxg to protocol://origin/d/dataset.cxg.
  TODO(cc) revisit special handling for local environment.
   */
  if (globals.API.local) {
    const dataRoot = "d";
    const deploymentId = bindDeploymentId(deploymentUrl);
    return `${window.location.origin}/${dataRoot}/${deploymentId}/`;
  }
  return deploymentUrl;
}

export function createExplorerUrl() {
  /*
  The current URL is passed as an "explorer URL" query string parameter to the dataset meta API. For environments where
  the current origin does not match the origin specified in the dataset deployment URLs (eg local, canary), update the
  origin to be the origin specified in the globals. Also update the data root for local environments.
  TODO(cc) revisit special handling for local and canary environments
   */
  const url = window.location.href;
  if (url.indexOf(globals.API.origin) === 0) {
    return url;
  }
  const { origin } = globals.API;
  if (globals.API.local) {
    const dataRoot = "e";
    const deploymentId = bindDeploymentId(url);
    return `${origin}${dataRoot}/${deploymentId}/`;
  }
  const dataRootAndDeploymentId = bindDataRootAndDeploymentId(url);
  return `${origin}${dataRootAndDeploymentId}/`;
}

export function sortDatasets(vm0, vm1) {
  /*
  Sort datasets by cell count, descending.
   */
  return (vm1.cell_count ?? 0) - (vm0.cell_count ?? 0);
}

function bindDataRootAndDeploymentId(pathOrUrl) {
  /*
  Read data root and deployment ID from path. Given "protocol://origin/x/any-alpha-numeric.cxg/", match on
  "/x/any-alpha-numeric.cxg/", group on "x/any-alpha-numeric.cxg".
   */
  const matches = pathOrUrl.match(/\/([a-z0-9_.-]+\/[a-z0-9_.-]+\.cxg)\//i);
  if (!matches || matches.length < 2) {
    if (globals.API.local) {
      return "e/03aa05fe-74dc-456f-ab2b-4b56d51f3ed0.cxg"; // TODO(cc) default data root and deployment ID for local (single mode)
    }
    // Expecting at least match and one capturing group
    throw new Error(
      `Unable to bind data root and deployment ID from "${pathOrUrl}"`
    );
  }
  return matches[1];
}

function bindDeploymentId(pathOrUrl) {
  /*
  read name of cxg from path. given "protocol://origin/x/any-alpha-numeric.cxg/", match on "/any-alpha-numeric.cxg/",
  group on "any-alpha-numeric.cxg".
   */
  const matches = pathOrUrl.match(/\/([a-z0-9_.-]*\.cxg)\//i);
  if (!matches || matches.length < 2) {
    if (globals.API.local) {
      return "03aa05fe-74dc-456f-ab2b-4b56d51f3ed0.cxg"; // TODO(cc) default dataset for local (single mode)
    }
    // Expecting at least match and one capturing group
    throw new Error(`Unable to bind deployment ID from "${pathOrUrl}"`);
  }
  return matches[1];
}
