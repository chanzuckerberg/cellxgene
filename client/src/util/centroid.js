import quantile from "./quantile";
import { memoize } from "./dataframe/util";

/*
  Centroid coordinate calculation
*/

/* Unused - please cleanup
const calcMeanCentroid = (world, annoName, annoValue, layoutDimNames) => {
  const centroid = { x: 0, y: 0, size: 0 };
  const annoArray = world.obsAnnotations.col(categoryName).asArray();
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

/* 
  calcMedianCentroid goes through a given metadata category
  fetches each cell's coordinates grouping by category value.

  It then calculates the median value and puts that in the array
*/

const calcMedianCentroid = (
  world,
  categoryName,
  layoutDimNames,
  categoricalSelection
) => {
  const categoryArray = world.obsAnnotations.col(categoryName).asArray();
  const layoutXArray = world.obsLayout.col(layoutDimNames[0]).asArray();
  const layoutYArray = world.obsLayout.col(layoutDimNames[1]).asArray();
  const coordinates = new Map();
  // Iterate over all the cells in the category
  for (let i = 0, len = categoryArray.length; i < len; i += 1) {
    const categoryValue = categoryArray[i];

    // Get the index of the categoryValue within the category
    const categoryValueIndex = categoricalSelection[
      categoryName
    ].categoryValueIndices.get(categoryValue);

    // Get the number of cells which are in the category value
    const numInCategoryValue =
      categoricalSelection[categoryName].categoryValueCounts[
        categoryValueIndex
      ];

    // Create/fetcht the valueArray,
    // which is what the key points to in the `coordinates` hashmap
    const valueArray = coordinates.get(categoryValue) || [
      false, // hasFinte
      0, // index
      new Float32Array(numInCategoryValue), // x coorindates
      new Float32Array(numInCategoryValue) // y coorindates
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

  // Iterate over the recently created map
  const iter = coordinates.entries();
  let pair = iter.next().value;
  let value;
  let key;
  while (pair) {
    key = pair[0];
    value = pair[1];
    // If there is a finite coordinate for the category value
    if (value[0]) {
      // Find the median x and y coordinate
      // and insert them into the first two indices
      value[0] = quantile([0.5], value[2])[0];
      value[1] = quantile([0.5], value[3])[0];
      // Remove the last two elements (where the arrays of coordinates were)
      value.pop();
      value.pop();
    } else {
      // remove the entry if not
      coordinates.delete(key);
    }
    pair = iter.next().value;
  }
  // return the map: categoricalValue -> [medianXCoordinate, medianYCoordinate]
  return coordinates;
};

// A simple function to hash the parameters
// (not 100% on world hash, Bruce will have to check this one out)
const hashMedianCentroid = (world, categoryName, layoutDimNames) => {
  return `${world.varAnnotations.__id}+${world.obsLayout.__id}+${
    world.varData.__id
  }::${categoryName}:${layoutDimNames}`;
};
// export the mmemoized calculation function
export default memoize(calcMedianCentroid, hashMedianCentroid);
