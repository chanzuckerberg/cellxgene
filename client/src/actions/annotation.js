/*
Action creators for user annotation
*/
import { unassignedCategoryLabel } from "../globals";
import { AnnoMatrixObsCrossfilter } from "../util/annoMatrix";

export const annotationCreateCategoryAction = (
  newCategoryName,
  categoryToDuplicate
) => async (dispatch, getState) => {
  /*
  Add a new user-created category to the obs annotations.
  Arguments:
    newCategoryName - string name for the category.
    categoryToDuplicate - obs category to use for initial values, or null.
  */
  const { annoMatrix: prevAnnoMatrix } = getState();
  const { schema } = prevAnnoMatrix;

  /* name must be a string,  non-zero length */
  if (typeof newCategoryName !== "string" || newCategoryName.length === 0)
    throw new Error("user annotations require string name");

  /* ensure the name isn't already in use! */
  if (schema.annotations.obsByName[newCategoryName])
    throw new Error("name collision on annotation category create");

  /* if we are duplicating a category, get it */
  let initialValue;
  let categories;
  if (categoryToDuplicate) {
    const catDupSchema = schema.annotations.obsByName[categoryToDuplicate];
    const catDupType = catDupSchema?.type;
    if (catDupType !== "string" && catDupType !== "categorical")
      throw new Error("categoryToDuplicate does not exist or has invalid type");

    const catToDupDf = await prevAnnoMatrix.fetch("obs", categoryToDuplicate);
    const col = catToDupDf.col(categoryToDuplicate);
    initialValue = col.asArray();
    ({ categories } = col.summarize());
  } else {
    initialValue = unassignedCategoryLabel;
    categories = [unassignedCategoryLabel];
  }

  const annoMatrix = prevAnnoMatrix.addObsColumn(
    {
      name: newCategoryName,
      type: "categorical",
      categories,
      writable: true,
    },
    Array,
    initialValue
  );

  const obsCrossfilter = new AnnoMatrixObsCrossfilter(annoMatrix);
  dispatch({
    type: "annotation: create category",
    data: newCategoryName,
    categoryToDuplicate,
    annoMatrix,
    obsCrossfilter,
  });
};

export const annotationRenameCategoryAction = (
  oldCategoryName,
  newCategoryName
) => (dispatch, getState) => {
  const { annoMatrix: prevAnnoMatrix } = getState();

  /* name must be a string,  non-zero length */
  if (typeof newCategoryName !== "string" || newCategoryName.length === 0)
    throw new Error("user annotations require string name");

  if (oldCategoryName === newCategoryName) return;

  const annoMatrix = prevAnnoMatrix.renameObsColumn(
    oldCategoryName,
    newCategoryName
  );
  const obsCrossfilter = new AnnoMatrixObsCrossfilter(annoMatrix);

  dispatch({
    type: "annotation: category edited",
    annoMatrix,
    obsCrossfilter,
    metadataField: oldCategoryName,
    newCategoryText: newCategoryName,
    data: newCategoryName,
  });
};

export const annotationDeleteCategoryAction = (categoryName) => (
  dispatch,
  getState
) => {
  const { annoMatrix: prevAnnoMatrix } = getState();
  const { schema } = prevAnnoMatrix;
  const colSchema = schema.annotations.obsByName[categoryName];

  if (!colSchema?.writable) throw new Error("unable to delete annotation");

  const annoMatrix = prevAnnoMatrix.dropObsColumn(categoryName);
  const obsCrossfilter = new AnnoMatrixObsCrossfilter(annoMatrix);

  dispatch({
    type: "annotation: delete category",
    annoMatrix,
    obsCrossfilter,
    metadataField: categoryName,
  });
};
