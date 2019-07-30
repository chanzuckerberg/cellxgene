import quantile from "./quantile";

/*
  Centroid coordinate calculation
*/

/* Unused - please cleanup
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
*/

const calcMedianCentroid = (world, annoName, layoutDimNames) => {
  const annoArray = world.obsAnnotations.col(annoName).asArray();
  const layoutXArray = world.obsLayout.col(layoutDimNames[0]).asArray();
  const layoutYArray = world.obsLayout.col(layoutDimNames[1]).asArray();

  const coordinates = new Map();

  for (let i = 0, len = annoArray.length; i < len; i += 1) {
    const valueArray = coordinates.get(annoArray[i]) || [false, [], []];
    let hasFinite = valueArray[0];
    hasFinite =
      Number.isFinite(layoutXArray[i]) || Number.isFinite(layoutYArray[i])
        ? true
        : hasFinite;
    valueArray[0] = hasFinite;
    valueArray[1].push(layoutXArray[i]);
    valueArray[2].push(layoutYArray[i]);
    coordinates.set(annoArray[i], valueArray);
  }

  const centroidCoordinates = [];

  const iter = coordinates.entries();

  let { value } = iter.next();
  while (value) {
    centroidCoordinates.push([
      value[0],
      ...quantile([0.5], new Float64Array(value[1][1])),
      ...quantile([0.5], new Float64Array(value[1][2]))
    ]);
    ({ value } = iter.next());
  }
  return centroidCoordinates;
};

const hashMedianCentroid = (world, annoName, layoutDimNames) => {
};

export default calcMedianCentroid;
