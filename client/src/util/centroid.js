import quantile from "./quantile";
import { memoize } from "./dataframe/util";
import { unassignedCategoryLabel } from "../globals";

/*
  Centroid coordinate calculation
*/

/* 
  calcMedianCentroid goes through a given metadata category
  fetches each cell's coordinates grouping by category value.

  It then calculates the median value and puts that in the array
*/

const calcMedianCentroid = (
  obsAnnotations,
  obsLayout,
  categoryName,
  layoutDimNames,
  categoricalSelection,
  schemaObsByName
) => {
  const categoryArray = obsAnnotations.col(categoryName).asArray();

  const layoutXArray = obsLayout.col(layoutDimNames[0]).asArray();
  const layoutYArray = obsLayout.col(layoutDimNames[1]).asArray();
  const coordinates = new Map();

  // Iterate over all the cells in the category
  for (let i = 0, len = categoryArray.length; i < len; i += 1) {
    const categoryValue = categoryArray[i];

    // Get the index of the categoryValue within the category
    // If the category is truncated and this value is removed,
    // it will not be assigned a category value and will not be
    // labeled on the graph
    const categoryValueIndex = categoricalSelection[
      categoryName
    ].categoryValueIndices.get(categoryValue);

    // Check to see if the current category is a user created annotation
    // if the user created this category, do not create a label for the `unassigned` value
    const isUserAnno = schemaObsByName[categoryName].writable;

    if (
      categoryValueIndex !== undefined &&
      !(isUserAnno && categoryValue === unassignedCategoryLabel)
    ) {
      // Get the number of cells which are in the category value
      const numInCategoryValue =
        categoricalSelection[categoryName].categoryValueCounts[
          categoryValueIndex
        ];

      // Create/fetch the valueArray,
      // which is what the key points to in the `coordinates` hashmap
      const valueArray = coordinates.get(categoryValue) || [
        false, // hasFinite
        0, // index
        new Float32Array(numInCategoryValue), // x coordinates
        new Float32Array(numInCategoryValue) // y coordinates
      ];
      const index = valueArray[1];
      let hasFinite = valueArray[0];

      hasFinite =
        Number.isFinite(layoutXArray[i]) && Number.isFinite(layoutYArray[i])
          ? true
          : hasFinite;
      valueArray[0] = hasFinite;
      valueArray[1] = index + 1;
      valueArray[2][index] = layoutXArray[i];
      valueArray[3][index] = layoutYArray[i];

      coordinates.set(categoryValue, valueArray);
    }
  }

  // Iterate over the recently created map
  coordinates.forEach((value, key) => {
    // If there are coordinates for this cateogrical value,
    // and there is a finite coordinate for the category value
    if (value[2].length > 0 && value[3].length > 0 && value[0]) {
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
  });
  // return the map: categoricalValue -> [medianXCoordinate, medianYCoordinate]
  return coordinates;
};

// A simple function to hash the parameters
// (not 100% on world hash, Bruce will have to check this one out)
const hashMedianCentroid = (
  obsAnnotations,
  obsLayout,
  categoryName,
  layoutDimNames,
  categorySelection,
  schemaObsByName
) => {
  return `${obsAnnotations.__id}+${
    obsLayout.__id
  }:${categoryName}:${layoutDimNames}:${Object.keys(
    categorySelection
  )}:${Object.keys(schemaObsByName)}`;
};
// export the mmemoized calculation function
export default memoize(calcMedianCentroid, hashMedianCentroid);
