/*
Helper functions for querying and binding Portal collections.
*/

export function lookupDatasetIdByPath(path, collections) {
  /*
  determine current dataset ID by mapping from URL path parameters to dataset values in specified set of collections 
  TODO(cc) temp workaround until dataset ID is added to response
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

function bindDeploymentId(path) {
  /*
  read name of cxg from path. given "protocol://hostname/x/any-alpha-numeric.cxg/", match on "/any-alpha-numeric.cxg/".
   */
  const matches = path.match(/\/([a-z0-9_.-]*\.cxg)\//i);
  if (!matches || matches.length < 2) {
    // Expecting at least match and one capturing group
    throw new Error("Unable to bind deployment ID from URL");
  }
  return matches[1];
}

function findDatasetIdByDeploymentId(deploymentId, collections) {
  /**
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
