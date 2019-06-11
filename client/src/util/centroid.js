/*
  Centroid coordinate calculation


*/


const calcMeanCentroid = (world, annoName, annoValue, layoutDimNames) => {
  const centroid = [0, 0, 0];
  const annoArray = world.obsAnnotations.col(annoName).asArray();
  const layoutXArray = world.obsLayout.col(layoutDimNames[0]).asArray();
  const layoutYArray = world.obsLayout.col(layoutDimNames[1]).asArray();

  for (let i = 0, len = annoArray.length; i < len; i += 1) {
    if (annoArray[i] === annoValue) {
      centroid[0] += layoutXArray[i];
      centroid[1] += layoutYArray[i];
      centroid[2] += 1;
    }
  }

  if (centroid[2] !== 0) {
    centroid[0] /= centroid[2];
    centroid[1] /= centroid[2];
    centroid[2] /= world.nObs;
  }

  return centroid;
};

export default calcMeanCentroid;
