import quantile from "./quantile";
import { memoize } from "./dataframe/util";

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

const calcMedianCentroid = (
  world,
  annoName,
  layoutDimNames,
  categoricalSelection
) => {
  const annoArray = world.obsAnnotations.col(annoName).asArray();
  const layoutXArray = world.obsLayout.col(layoutDimNames[0]).asArray();
  const layoutYArray = world.obsLayout.col(layoutDimNames[1]).asArray();

  const coordinates = new Map(); // map: annoValue -> [hasFinite, counter, [x Float32], [y Float32]]
  for (let i = 0, len = annoArray.length; i < len; i += 1) {
    const categoryValue = annoArray[i];

    const categoryValueIndex = categoricalSelection[
      annoName
    ].categoryValueIndices.get(categoryValue);

    const numInAnno =
      categoricalSelection[annoName].categoryValueCounts[categoryValueIndex];

    const valueArray = coordinates.get(categoryValue) || [
      false,
      0,
      new Float32Array(numInAnno),
      new Float32Array(numInAnno)
    ];
    const index = valueArray[1];
    let hasFinite = valueArray[0];
    hasFinite =
      Number.isFinite(layoutXArray[i]) || Number.isFinite(layoutYArray[i])
        ? true
        : hasFinite;
    valueArray[0] = hasFinite;
    valueArray[1] = index + 1;
    valueArray[2][index] = layoutXArray[i];
    valueArray[3][index] = layoutYArray[i];
    coordinates.set(categoryValue, valueArray);
  }

  const iter = coordinates.entries();
  let pair = iter.next().value;
  let value;
  let key;
  while (pair) {
    key = pair[0];
    value = pair[1];
    if (value[0]) {
      value[0] = quantile([0.5], value[2])[0];
      value[1] = quantile([0.5], value[3])[0];
      value.pop();
      value.pop();
    } else {
      coordinates.delete(key);
    }
    pair = iter.next().value;
  }

  return coordinates;
};

const hashMedianCentroid = (world, annoName, layoutDimNames) => {
  return `${world.varAnnotations.__id}+${world.obsLayout.__id}+${
    world.varData.__id
  }::${annoName}:${layoutDimNames}`;
};
export default memoize(calcMedianCentroid, hashMedianCentroid);
