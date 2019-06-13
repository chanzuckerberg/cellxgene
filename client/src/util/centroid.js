import quantile from "./quantile";

/*
  Centroid coordinate calculation
*/
const calcMeanCentroid = (world, annoName, annoValue, layoutDimNames) => {
  const centroid = { x: 0, y: 0, size: 0 };
  const annoArray = world.obsAnnotations.col(annoName).asArray();
  const layoutXArray = world.obsLayout.col(layoutDimNames[0]).asArray();
  const layoutYArray = world.obsLayout.col(layoutDimNames[1]).asArray();

  for (let i = 0, len = annoArray.length; i < len; i += 1) {
    if (annoArray[i] === annoValue) {
      centroid.x += layoutXArray[i];
      centroid.y += layoutYArray[i];
      centroid.size += 1;
    }
  }

  if (centroid[2] !== 0) {
    centroid.x /= centroid.size;
    centroid.y /= centroid.size;
  }

  return [centroid.x, centroid.y];
};

const calcMedianCentroid = (world, annoName, annoValue, layoutDimNames) => {
  const centroidX = [];
  const centroidY = [];
  let hasFinite = false;

  const annoArray = world.obsAnnotations.col(annoName).asArray();
  const layoutXArray = world.obsLayout.col(layoutDimNames[0]).asArray();
  const layoutYArray = world.obsLayout.col(layoutDimNames[1]).asArray();

  for (let i = 0, len = annoArray.length; i < len; i += 1) {
    if (annoArray[i] === annoValue) {
      hasFinite = Number.isFinite(annoArray[i]) ? true : hasFinite;
      centroidX.push(layoutXArray[i]);
      centroidY.push(layoutYArray[i]);
    }
  }
  if (hasFinite) {
    const medianX = quantile([0.5], Float64Array.from(centroidX));
    const medianY = quantile([0.5], Float64Array.from(centroidY));

    return [medianX, medianY];
  }

  return null;
};

export default calcMedianCentroid;
