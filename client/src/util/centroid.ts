import quantile from "./quantile";
import { memoize } from "./dataframe/util";
import { Dataframe } from "./dataframe";
import { unassignedCategoryLabel } from "../globals";
import {
  createCategorySummaryFromDfCol,
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
const getCoordinatesByLabel = (
  // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
  schema: any,
  categoryName: string,
  categoryDf: Dataframe,
  // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
  layoutChoice: any,
  layoutDf: Dataframe
) => {
  const coordsByCategoryLabel = new Map();
  // If the coloredBy is not a categorical col
  if (!isSelectableCategoryName(schema, categoryName)) {
    return coordsByCategoryLabel;
  }

  const categoryArray = categoryDf.col(categoryName).asArray();
  const layoutDimNames = layoutChoice.currentDimNames;
  const layoutXArray = layoutDf.col(layoutDimNames[0]).asArray();
  const layoutYArray = layoutDf.col(layoutDimNames[1]).asArray();

  const categorySummary = createCategorySummaryFromDfCol(
    categoryDf.col(categoryName),
    schema.annotations.obsByName[categoryName]
  );

  const { isUserAnno, categoryValueIndices, categoryValueCounts } =
    categorySummary;

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
        // @ts-expect-error ts-migrate(2538) FIXME: Blocked by StateManager/ControlsHelpers
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

const calcMedianCentroid = (
  // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
  schema: any,
  categoryName: string,
  categoryDf: Dataframe,
  // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
  layoutChoice: any,
  layoutDf: Dataframe
) => {
  // generate a map describing the coordinates for each label within the given category
  const dataMap = getCoordinatesByLabel(
    schema,
    categoryName,
    categoryDf,
    layoutChoice,
    layoutDf
  );

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
const hashMedianCentroid = (
  // @ts-expect-error ts-migrate(6133) FIXME: 'schema' is declared but its value is never read.
  // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
  schema: any,
  categoryName: string,
  categoryDf: Dataframe,
  // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
  layoutChoice: any,
  layoutDf: Dataframe
): string => {
  const category = categoryDf.col(categoryName);
  const layoutDimNames = layoutChoice.currentDimNames;
  const layoutX = layoutDf.col(layoutDimNames[0]);
  const layoutY = layoutDf.col(layoutDimNames[1]);
  return `${category.__id}+${layoutX.__id}:${layoutY.__id}`;
};
// export the memoized calculation function
export default memoize(calcMedianCentroid, hashMedianCentroid);
