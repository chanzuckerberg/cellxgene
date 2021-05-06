/*
Helper functions for querying and binding Portal collections.
*/

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
