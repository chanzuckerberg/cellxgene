/*
Action creators for user annotation
*/
import _ from "lodash";
import pako from "pako";
import * as globals from "../globals";
import { MatrixFBS, AnnotationsHelpers } from "../util/stateManager";

const { isUserAnnotation } = AnnotationsHelpers;

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
  let newSchema;
  let ctor;
  if (categoryToDuplicate) {
    /* if we are duplicating a category, retrieve it */
    const catDupSchema = schema.annotations.obsByName[categoryToDuplicate];
    const catDupType = catDupSchema?.type;
    if (catDupType !== "string" && catDupType !== "categorical")
      throw new Error("categoryToDuplicate does not exist or has invalid type");

    const catToDupDf = await prevAnnoMatrix
      .base()
      .fetch("obs", categoryToDuplicate);
    const col = catToDupDf.col(categoryToDuplicate);
    initialValue = col.asArray();
    const { categories } = col.summarizeCategorical();
    // all user-created annotations must have the unassigned category
    if (!categories.includes(globals.unassignedCategoryLabel)) {
      categories.push(globals.unassignedCategoryLabel);
    }
    ctor = initialValue.constructor;
    newSchema = {
      ...catDupSchema,
      name: newCategoryName,
      categories,
      writable: true,
    };
  } else {
    /* else assign to the standard default value */
    initialValue = globals.unassignedCategoryLabel;
    ctor = Array;
    newSchema = {
      name: newCategoryName,
      type: "categorical",
      categories: [globals.unassignedCategoryLabel],
      writable: true,
    };
  }

  const obsCrossfilter = prevObsCrossfilter.addObsColumn(
    newSchema,
    ctor,
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
) => async (dispatch, getState) => {
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
    obsCrossfilter = await obsCrossfilter.setObsColumnValues(
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

export const annotationDeleteLabelFromCategory = (
  categoryName,
  labelName
) => async (dispatch, getState) => {
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

  const obsCrossfilter = await prevObsCrossfilter.removeObsAnnoCategory(
    categoryName,
    labelName,
    globals.unassignedCategoryLabel
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
) => async (dispatch, getState) => {
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

  let obsCrossfilter = await prevObsCrossfilter.resetObsColumnValues(
    categoryName,
    oldLabelName,
    newLabelName
  );
  obsCrossfilter = await obsCrossfilter.removeObsAnnoCategory(
    categoryName,
    oldLabelName,
    globals.unassignedCategoryLabel
  );

  dispatch({
    type: "annotation: label edited",
    editedLabel: newLabelName,
    metadataField: categoryName,
    label: oldLabelName,
    annoMatrix: obsCrossfilter.annoMatrix,
    obsCrossfilter,
  });
};

export const annotationLabelCurrentSelection = (
  categoryName,
  labelName
) => async (dispatch, getState) => {
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

  const obsCrossfilter = await prevObsCrossfilter.setObsColumnValues(
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

function writableAnnotations(annoMatrix) {
  return annoMatrix.schema.annotations.obs.columns
    .filter((s) => s.writable)
    .map((s) => s.name);
}

export const needToSaveObsAnnotations = (annoMatrix, lastSavedAnnoMatrix) => {
  /*
  Return true if there are LIKELY user-defined annotation modifications between the two
  annoMatrices.  Technically not an action creator, but intimately intertwined
  with the save process.

  Two conditions will trigger a need to save:
    * the collection of user-defined columns have changed
    * the contents of the user-defined columns have change
  */

  annoMatrix = annoMatrix.base();

  // if the annoMatrix hasn't changed, we are guaranteed no changes to the matrix schema or contents.
  if (annoMatrix === lastSavedAnnoMatrix) return false;

  // if the schema has changed, we need to save
  const currentWritable = writableAnnotations(annoMatrix);
  if (_.difference(currentWritable, writableAnnotations(lastSavedAnnoMatrix))) {
    return true;
  }

  // no schema changes; check for change in contents
  return currentWritable.some(
    (col) => annoMatrix.col(col) !== lastSavedAnnoMatrix.col(col)
  );
};

export const saveObsAnnotationsAction = () => async (dispatch, getState) => {
  /*
  Save the user-created obs annotations IF any have changed.
  */
  const state = getState();
  const { annotations, autosave } = state;
  const { dataCollectionNameIsReadOnly, dataCollectionName } = annotations;
  const { lastSavedAnnoMatrix, saveInProgress } = autosave;

  const annoMatrix = state.annoMatrix.base();

  if (saveInProgress || annoMatrix === lastSavedAnnoMatrix) return;
  if (!needToSaveObsAnnotations(annoMatrix, lastSavedAnnoMatrix)) {
    dispatch({
      type: "writable obs annotations - save complete",
      lastSavedAnnoMatrix: annoMatrix,
    });
    return;
  }

  /*
  Else, we really do need to save
  */

  dispatch({
    type: "writable obs annotations - save started",
  });

  const df = await annoMatrix.fetch("obs", writableAnnotations(annoMatrix));
  const matrix = MatrixFBS.encodeMatrixFBS(df);
  const compressedMatrix = pako.deflate(matrix);
  try {
    const queryString =
      !dataCollectionNameIsReadOnly && !!dataCollectionName
        ? `?annotation-collection-name=${encodeURIComponent(
            dataCollectionName
          )}`
        : "";
    const res = await fetch(
      `${globals.API.prefix}${globals.API.version}annotations/obs${queryString}`,
      {
        method: "PUT",
        body: compressedMatrix,
        headers: new Headers({
          "Content-Type": "application/octet-stream",
        }),
        credentials: "include",
      }
    );
    if (res.ok) {
      dispatch({
        type: "writable obs annotations - save complete",
        lastSavedAnnoMatrix: annoMatrix,
      });
    } else {
      dispatch({
        type: "writable obs annotations - save error",
        message: `HTTP error ${res.status} - ${res.statusText}`,
        res,
      });
    }
  } catch (error) {
    dispatch({
      type: "writable obs annotations - save error",
      message: error.toString(),
      error,
    });
  }
};

export const saveGenesetsAction = () => async (dispatch, getState) => {
  const state = getState();

  // bail if gene sets not available, or in readonly mode.
  const { config } = state;
  const { lastTid, genesets } = state.genesets;

  const genesetsAreAvailable =
    config?.parameters?.annotations_genesets ?? false;
  const genesetsReadonly =
    config?.parameters?.annotations_genesets_readonly ?? true;
  if (!genesetsAreAvailable || genesetsReadonly) {
    // our non-save was completed!
    return dispatch({
      type: "autosave: genesets complete",
      lastSavedGenesets: genesets,
    });
  }

  dispatch({
    type: "autosave: genesets started",
  });

  /* Create the JSON OTA data structure */
  const tid = (lastTid ?? 0) + 1;
  const genesetsOTA = [];
  for (const [name, gs] of genesets) {
    const genes = [];
    for (const g of gs.genes.values()) {
      genes.push({
        gene_symbol: g.geneSymbol,
        gene_description: g.geneDescription,
      });
    }
    genesetsOTA.push({
      geneset_name: name,
      geneset_description: gs.genesetDescription,
      genes,
    });
  }
  const ota = {
    tid,
    genesets: genesetsOTA,
  };

  /* Save to server */
  try {
    const {
      dataCollectionNameIsReadOnly,
      dataCollectionName,
    } = state.annotations;
    const queryString =
      !dataCollectionNameIsReadOnly && !!dataCollectionName
        ? `?annotation-collection-name=${encodeURIComponent(
            dataCollectionName
          )}`
        : "";

    const res = await fetch(
      `${globals.API.prefix}${globals.API.version}genesets${queryString}`,
      {
        method: "PUT",
        headers: new Headers({
          Accept: "application/json",
          "Content-Type": "application/json",
        }),
        body: JSON.stringify(ota),
        credentials: "include",
      }
    );
    if (!res.ok) {
      return dispatch({
        type: "autosave: genesets error",
        message: `HTTP error ${res.status} - ${res.statusText}`,
        res,
      });
    }
    return Promise.all([
      dispatch({
        type: "autosave: genesets complete",
        lastSavedGenesets: genesets,
      }),
      dispatch({
        type: "geneset: set tid",
        tid,
      }),
    ]);
  } catch (error) {
    return dispatch({
      type: "autosave: genesets error",
      message: error.toString(),
      error,
    });
  }
};
