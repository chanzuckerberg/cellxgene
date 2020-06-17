/*
Action creators for user annotation
*/
import { unassignedCategoryLabel } from "../globals";

function isUserAnnotation(annoMatrix, obsCatName) {
  // true if writable annotation; false if not or doesn't exist.
  return annoMatrix.schema.annotations.obsByName[obsCatName]?.writable ?? false;
}

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
  const {
    annoMatrix: prevAnnoMatrix,
    obsCrossfilter: prevObsCrossfilter,
  } = getState();
  if (!prevAnnoMatrix || !prevObsCrossfilter) return;
  const { schema } = prevAnnoMatrix;

  /* name must be a string,  non-zero length */
  if (typeof newCategoryName !== "string" || newCategoryName.length === 0)
    throw new Error("user annotations require string name");

  /* ensure the name isn't already in use! */
  if (schema.annotations.obsByName[newCategoryName])
    throw new Error("name collision on annotation category create");

  let initialValue;
  let categories;
  if (categoryToDuplicate) {
    /* if we are duplicating a category, retrieve it */
    const catDupSchema = schema.annotations.obsByName[categoryToDuplicate];
    const catDupType = catDupSchema?.type;
    if (catDupType !== "string" && catDupType !== "categorical")
      throw new Error("categoryToDuplicate does not exist or has invalid type");

    const catToDupDf = await prevAnnoMatrix.fetch("obs", categoryToDuplicate);
    const col = catToDupDf.col(categoryToDuplicate);
    initialValue = col.asArray();
    ({ categories } = col.summarize());
  } else {
    /* else assign to the standard default value */
    initialValue = unassignedCategoryLabel;
    categories = [unassignedCategoryLabel];
  }

  const obsCrossfilter = prevObsCrossfilter.addObsColumn(
    {
      name: newCategoryName,
      type: "categorical",
      categories,
      writable: true,
    },
    Array,
    initialValue
  );

  dispatch({
    type: "annotation: create category",
    data: newCategoryName,
    categoryToDuplicate,
    annoMatrix: obsCrossfilter.annoMatrix,
    obsCrossfilter,
  });
};

export const annotationRenameCategoryAction = (
  oldCategoryName,
  newCategoryName
) => (dispatch, getState) => {
  /*
  Rename a user-created annotation category
  */
  const {
    annoMatrix: prevAnnoMatrix,
    obsCrossfilter: prevObsCrossfilter,
  } = getState();
  if (!prevAnnoMatrix || !prevObsCrossfilter) return;
  if (!isUserAnnotation(prevAnnoMatrix, oldCategoryName))
    throw new Error("not a user annotation");

  /* name must be a string,  non-zero length */
  if (typeof newCategoryName !== "string" || newCategoryName.length === 0)
    throw new Error("user annotations require string name");

  if (oldCategoryName === newCategoryName) return;

  const obsCrossfilter = prevObsCrossfilter.renameObsColumn(
    oldCategoryName,
    newCategoryName
  );

  dispatch({
    type: "annotation: category edited",
    annoMatrix: obsCrossfilter.annoMatrix,
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
  /*
  Delete a user-created category
  */
  const {
    annoMatrix: prevAnnoMatrix,
    obsCrossfilter: prevObsCrossfilter,
  } = getState();
  if (!prevAnnoMatrix || !prevObsCrossfilter) return;
  if (!isUserAnnotation(prevAnnoMatrix, categoryName))
    throw new Error("not a user annotation");

  const obsCrossfilter = prevObsCrossfilter.dropObsColumn(categoryName);
  dispatch({
    type: "annotation: delete category",
    annoMatrix: obsCrossfilter.annoMatrix,
    obsCrossfilter,
    metadataField: categoryName,
  });
};

export const annotationCreateLabelInCategory = (
  categoryName,
  labelName,
  assignSelected
) => (dispatch, getState) => {
  /*
  Add a new label to a user-defined category.  If assignSelected is true, assign 
  the label to all currently selected cells.
  */
  const {
    annoMatrix: prevAnnoMatrix,
    obsCrossfilter: prevObsCrossfilter,
  } = getState();
  if (!prevAnnoMatrix || !prevObsCrossfilter) return;
  if (!isUserAnnotation(prevAnnoMatrix, categoryName))
    throw new Error("not a user annotation");

  let obsCrossfilter = prevObsCrossfilter.addObsAnnoCategory(
    categoryName,
    labelName
  );
  if (assignSelected) {
    obsCrossfilter = obsCrossfilter.setObsColumnValues(
      categoryName,
      prevObsCrossfilter.allSelectedLabels(),
      labelName
    );
  }

  dispatch({
    type: "annotation: add new label to category",
    annoMatrix: obsCrossfilter.annoMatrix,
    obsCrossfilter,
    metadataField: categoryName,
    newLabelText: labelName,
    assignSelectedCells: assignSelected,
  });
};

export const annotationDeleteLabelFromCategory = (categoryName, labelName) => (
  dispatch,
  getState
) => {
  /*
  delete a label from a user-defined category
  */
  const {
    annoMatrix: prevAnnoMatrix,
    obsCrossfilter: prevObsCrossfilter,
  } = getState();
  if (!prevAnnoMatrix || !prevObsCrossfilter) return;
  if (!isUserAnnotation(prevAnnoMatrix, categoryName))
    throw new Error("not a user annotation");

  const obsCrossfilter = prevObsCrossfilter.removeObsAnnoCategory(
    categoryName,
    labelName,
    unassignedCategoryLabel
  );

  dispatch({
    type: "annotation: delete label",
    metadataField: categoryName,
    label: labelName,
    annoMatrix: obsCrossfilter.annoMatrix,
    obsCrossfilter,
  });
};

export const annotationRenameLabelInCategory = (
  categoryName,
  oldLabelName,
  newLabelName
) => (dispatch, getState) => {
  /*
  label name change
  */
  const {
    annoMatrix: prevAnnoMatrix,
    obsCrossfilter: prevObsCrossfilter,
  } = getState();
  if (!prevAnnoMatrix || !prevObsCrossfilter) return;
  if (!isUserAnnotation(prevAnnoMatrix, categoryName))
    throw new Error("not a user annotation");

  const obsCrossfilter = prevObsCrossfilter
    .resetObsColumnValues(categoryName, oldLabelName, newLabelName)
    .removeObsAnnoCategory(categoryName, oldLabelName, unassignedCategoryLabel);

  dispatch({
    type: "annotation: label edited",
    editedLabel: newLabelName,
    metadataField: categoryName,
    label: oldLabelName,
    annoMatrix: obsCrossfilter.annoMatrix,
    obsCrossfilter,
  });
};

export const annotationLabelCurrentSelection = (categoryName, labelName) => (
  dispatch,
  getState
) => {
  /*
  set the label on all currently selected
  */
  const {
    annoMatrix: prevAnnoMatrix,
    obsCrossfilter: prevObsCrossfilter,
  } = getState();
  if (!prevAnnoMatrix || !prevObsCrossfilter) return;
  if (!isUserAnnotation(prevAnnoMatrix, categoryName))
    throw new Error("not a user annotation");

  const obsCrossfilter = prevObsCrossfilter.setObsColumnValues(
    categoryName,
    prevObsCrossfilter.allSelectedLabels(),
    labelName
  );

  dispatch({
    type: "annotation: label current cell selection",
    metadataField: categoryName,
    label: labelName,
    obsCrossfilter,
    annoMatrix: obsCrossfilter.annoMatrix,
  });
};
