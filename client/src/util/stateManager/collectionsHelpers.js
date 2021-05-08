/*
Helper functions for querying and binding Portal collections.
*/

const LOCATION_ORIGIN = window.location.origin;

export function createCollectionsViewModel(collections) {
  /*
  update origin and data root of dataset deployment URLs to format expected for current environment, sort datasets
  by created date.
  TODO(cc) temp while we read collections from static JSON file
   */
  return collections.map((collection) => {
    return {
      ...collection,
      datasets: createDatasetsViewModel(collection.datasets),
    };
  });
}

export function lookupDatasetIdByPath(path, collections) {
  /*
  determine current dataset ID by mapping from URL path parameters to dataset values in specified set of collections .
  TODO(cc) temp until dataset ID is added to response
   */
  const deploymentId = bindDeploymentId(path);
  return findDatasetIdByDeploymentId(deploymentId, collections);
}

export function keyCollectionsByDatasetId(collections) {
  /*
  create convenience map facilitating select of collection for a given dataset.
  */
  return collections.reduce((collectionsByDatasetId, collection) => {
    (collection.datasets ?? []).forEach((dataset) => {
      collectionsByDatasetId.set(dataset.id, collection);
    });

    return collectionsByDatasetId;
  }, new Map());
}

function bindDeploymentId(pathOrUrl) {
  /*
  read name of cxg from path. given "protocol://hostname/x/any-alpha-numeric.cxg/", match on "/any-alpha-numeric.cxg/",
  group on "any-alpha-numeric.cxg".
   */
  const matches = pathOrUrl.match(/\/([a-z0-9_.-]*\.cxg)\//i);
  if (!matches || matches.length < 2) {
    // Expecting at least match and one capturing group
    throw new Error(`Unable to bind deployment ID from "${pathOrUrl}"`);
  }
  return matches[1];
}

function switchUrlOriginToCurrent(url) {
  /*
  swap out origin in the given URL to the current origin.
  TODO(cc) generalize data root per environment, use /\/[a-z0-9_.-]*\/[a-z0-9_.-]*\.cxg\//i to match data root and deployment ID.
   */
  const deploymentId = bindDeploymentId(url);
  return `${LOCATION_ORIGIN}/d/${deploymentId}/`;
}

function findDatasetIdByDeploymentId(deploymentId, collections) {
  /*
   return the ID of the dataset that has a deployment URL containing the specified deployment ID. 
   */
  let selectedDatasetId;
  for (const collection of collections) {
    for (const dataset of collection.datasets) {
      if (dataset.url.indexOf(deploymentId) >= 0) {
        selectedDatasetId = dataset.id;
        break;
      }
    }
    if (selectedDatasetId) {
      break;
    }
  }
  if (!selectedDatasetId) {
    throw new Error(
      `Unable to lookup dataset ID for deployment "${deploymentId}"`
    );
  }
  return selectedDatasetId;
}

function createDatasetsViewModel(datasets) {
  /*
   update origin of deployment URL of each dataset to be the current origin. sort by created date descending.
   */
  const viewModels = datasets.map((dataset) => {
    return {
      ...dataset,
      url: switchUrlOriginToCurrent(dataset.url),
    };
  });
  viewModels.sort((vm0, vm1) => {
    return vm0.createdAt - vm1.createdAt;
  });
  return viewModels;
}
