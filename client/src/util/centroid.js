import quantile from "./quantile";
import { memoize } from "./dataframe/util";
import { unassignedCategoryLabel } from "../globals";

/*
  Centroid coordinate calculation
  Calculates centroids for displaying
  In the case that a category is truncated, the truncated labels will not have
    centroids calculated
*/

/*
Generates a mapping of categorical values to data needed to calculate centroids
categoricalValue -> {
 length: int,
 holdsFinite: Boolean,
 xCoordinates: Float32Array,
 yCoordinates: Float32Array
}
*/
const getCoordinatesByCategoricalValues = (
  obsAnnotations,
  obsLayout,
  categoryName,
  layoutDimNames,
  categoricalSelection,
  schemaObsByName
) => {
  const coordsByCategoryLabel = new Map();

  const categoryArray = obsAnnotations.col(categoryName).asArray();

  const layoutXArray = obsLayout.col(layoutDimNames[0]).asArray();
  const layoutYArray = obsLayout.col(layoutDimNames[1]).asArray();

  const { categoryValueIndices, categoryValueCounts } =
    categoricalSelection?.[categoryName] || {};

  // If the coloredBy is not a categorical col
  if (categoryValueIndices === undefined) {
    return coordsByCategoryLabel;
  }

  // Check to see if the current category is a user created annotation
  const isUserAnno = schemaObsByName[categoryName].writable;

  // Iterate over all cells
  for (let i = 0, len = categoryArray.length; i < len; i += 1) {
    // Fetch the categorical value of the current cell
    const categoryValue = categoryArray[i];

    // Get the index of the categoryValue within the category
    const categoryValueIndex = categoryValueIndices.get(categoryValue);

    // If the category is truncated and this value is removed,
    //  it will not be assigned a category value and will not be
    //  labeled on the graph
    // If the user created this category,
    //  do not create a label for the `unassigned` value
    if (
      categoryValueIndex !== undefined &&
      !(isUserAnno && categoryValue === unassignedCategoryLabel)
    ) {
      // Create/fetch the scratchpad value
      let coords = coordsByCategoryLabel.get(categoryValue);
      if (coords === undefined) {
        // Get the number of cells which are in the categorical value
        const numInCategoricalValue = categoryValueCounts[categoryValueIndex];
        coords = {
          hasFinite: false,
          xCoordinates: new Float32Array(numInCategoricalValue),
          yCoordinates: new Float32Array(numInCategoricalValue),
          length: 0
        };
        coordsByCategoryLabel.set(categoryValue, coords);
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
  calcMedianCentroid calculates the median coordinates for categorical values in a given metadata field 

  categoricalValue -> [x-Coordinate, y-Coordinate]
*/

const calcMedianCentroid = (
  obsAnnotations,
  obsLayout,
  categoryName,
  layoutDimNames,
  categoricalSelection,
  schemaObsByName
) => {
  // generate a map describing the coordinates for each value within the given category
  const dataMap = getCoordinatesByCategoricalValues(
    obsAnnotations,
    obsLayout,
    categoryName,
    layoutDimNames,
    categoricalSelection,
    schemaObsByName
  );

  // categoricalValue => [medianXCoordinate, medianYCoordinate]
  const coordinates = new Map();

  // Iterate over the recently created map
  dataMap.forEach((value, key) => {
    // If there are coordinates for this categorical value,
    // and there is a finite coordinate for the category value
    if (value.length > 0 && value.hasFinite) {
      const calculatedCoordinates = [];

      // Find and store the median x and y coordinate
      calculatedCoordinates[0] = quantile([0.5], value.xCoordinates)[0];
      calculatedCoordinates[1] = quantile([0.5], value.yCoordinates)[0];

      coordinates.set(key, calculatedCoordinates);
    }
  });

  // return the map: categoricalValue -> [medianXCoordinate, medianYCoordinate]
  return coordinates;
};

// A simple function to hash the parameters
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
// export the memoized calculation function
export default memoize(calcMedianCentroid, hashMedianCentroid);
