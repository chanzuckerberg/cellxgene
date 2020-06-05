import quantile from "./quantile";
import { memoize } from "./dataframe/util";
import { unassignedCategoryLabel } from "../globals";
import {
  createCategorySummary,
  isSelectableCategoryName,
} from "./stateManager/controlsHelpers";

/*
  Centroid coordinate calculation
  Calculates centroids for displaying
  In the case that a category is truncated, the truncated labels will not have
    centroids calculated
*/

/*
Generates a mapping of labels to data needed to calculate centroids
label -> {
 length: int,
 holdsFinite: Boolean,
 xCoordinates: Float32Array,
 yCoordinates: Float32Array
}
*/
const getCoordinatesByLabel = (world, categoryName, layoutDimNames) => {
  const coordsByCategoryLabel = new Map();
  // If the coloredBy is not a categorical col
  if (!isSelectableCategoryName(world.schema, categoryName)) {
    return coordsByCategoryLabel;
  }

  const { obsAnnotations, obsLayout } = world;
  const categoryArray = obsAnnotations.col(categoryName).asArray();
  const layoutXArray = obsLayout.col(layoutDimNames[0]).asArray();
  const layoutYArray = obsLayout.col(layoutDimNames[1]).asArray();

  const categorySummary = createCategorySummary(world, categoryName);
  const {
    isUserAnno,
    categoryValueIndices,
    categoryValueCounts,
  } = categorySummary;

  // Iterate over all cells
  for (let i = 0, len = categoryArray.length; i < len; i += 1) {
    // Fetch the label of the current cell
    const label = categoryArray[i];

    // Get the index of the label within the category
    const labelIndex = categoryValueIndices.get(label);

    // If the category's labels are truncated and this label is removed,
    //  it will not be assigned a label and will not be
    //  labeled on the graph
    // If the user created this category,
    //  do not create a coord for the `unassigned` label
    if (
      labelIndex !== undefined &&
      !(isUserAnno && label === unassignedCategoryLabel)
    ) {
      // Create/fetch the scratchpad value
      let coords = coordsByCategoryLabel.get(label);
      if (coords === undefined) {
        // Get the number of cells which are in the label
        const numInLabel = categoryValueCounts[labelIndex];
        coords = {
          hasFinite: false,
          xCoordinates: new Float32Array(numInLabel),
          yCoordinates: new Float32Array(numInLabel),
          length: 0,
        };
        coordsByCategoryLabel.set(label, coords);
      }

      coords.hasFinite =
        coords.hasFinite ||
        (Number.isFinite(layoutXArray[i]) && Number.isFinite(layoutYArray[i]));

      const coordinatesLength = coords.length;

      coords.xCoordinates[coordinatesLength] = layoutXArray[i];
      coords.yCoordinates[coordinatesLength] = layoutYArray[i];

      coords.length = coordinatesLength + 1;
    }
  }
  return coordsByCategoryLabel;
};

/* 
  calcMedianCentroid calculates the median coordinates for labels in a given category

  label -> [x-Coordinate, y-Coordinate]
*/

const calcMedianCentroid = (world, categoryName, layoutDimNames) => {
  // generate a map describing the coordinates for each label within the given category
  const dataMap = getCoordinatesByLabel(world, categoryName, layoutDimNames);

  // label => [medianXCoordinate, medianYCoordinate]
  const coordinates = new Map();

  // Iterate over the recently created map
  dataMap.forEach((coords, label) => {
    // If there are coordinates for this label,
    // and there is a finite coordinate for the label
    if (coords.length > 0 && coords.hasFinite) {
      const calculatedCoordinates = [];

      // Find and store the median x and y coordinate
      calculatedCoordinates[0] = quantile([0.5], coords.xCoordinates)[0];
      calculatedCoordinates[1] = quantile([0.5], coords.yCoordinates)[0];

      coordinates.set(label, calculatedCoordinates);
    }
  });

  // return the map: label -> [medianXCoordinate, medianYCoordinate]
  return coordinates;
};

// A simple function to hash the parameters
const hashMedianCentroid = (world, categoryName, layoutDimNames) => {
  const { obsAnnotations, obsLayout } = world;
  return `${obsAnnotations.__id}+${obsLayout.__id}:${categoryName}:${layoutDimNames}`;
};
// export the memoized calculation function
export default memoize(calcMedianCentroid, hashMedianCentroid);
