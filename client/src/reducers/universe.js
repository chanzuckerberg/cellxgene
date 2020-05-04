import { unassignedCategoryLabel } from "../globals";
import {
  addObsAnnotations,
  addVarAnnotations,
  addObsLayout,
} from "../util/stateManager/universe";
import {
  World,
  ControlsHelpers,
  AnnotationsHelpers,
} from "../util/stateManager";

const Universe = (state = null, action, nextSharedState, prevSharedState) => {
  switch (action.type) {
    case "universe exists, but loading is still in progress": {
      const { universe } = action;
      return universe;
    }

    case "universe: column load success": {
      const { dim, dataframe } = action;
      switch (dim) {
        case "obsAnnotations": {
          return {
            ...state,
            ...addObsAnnotations(state, dataframe),
          };
        }
        case "varAnnotations": {
          return {
            ...state,
            ...addVarAnnotations(state, dataframe),
          };
        }
        case "obsLayout": {
          return {
            ...state,
            ...addObsLayout(state, dataframe),
          };
        }
        default: {
          throw new Error("action handler not implemented");
        }
      }
    }

    case "expression load success": {
      let { varData } = state;

      // Lazy load new expression data into the varData dataframe, if
      // not already present.
      //
      Object.entries(action.expressionData).forEach(([key, val]) => {
        // If not already in universe.varData, save entire expression column
        if (!varData.hasCol(key)) {
          varData = varData.withCol(key, val);
        }
      });

      // Prune size of varData "cache" if getting out of hand....
      //
      const { userDefinedGenes, diffexpGenes } = prevSharedState;
      const allTheGenesWeNeed = [
        ...new Set(
          [userDefinedGenes, diffexpGenes, Object.keys(action.expressionData)]
            .filter((ele) => ele)
            .flat()
        ),
      ];
      varData = ControlsHelpers.pruneVarDataCache(varData, allTheGenesWeNeed);
      return { ...state, varData };
    }

    case "annotation: create category": {
      /* create a new annotation category, with all values set to 'unassigned' */
      const name = action.data;
      const { categoryToDuplicate } = action;

      /* name must be a string,  non-zero length */
      if (typeof name !== "string" || name.length === 0)
        throw new Error("user annotations require string name");

      /* ensure the name isn't already in use! */
      if (state.obsAnnotations.hasCol(name))
        throw new Error("name collision on annotation category create");

      /* ensure the duplicate col exists */
      if (
        categoryToDuplicate &&
        !state.obsAnnotations.hasCol(categoryToDuplicate)
      )
        throw new Error("categoryToDuplicate does not exist");

      let schema;
      let data;
      if (categoryToDuplicate) {
        /* duplicate the named annotation */
        schema = AnnotationsHelpers.dupObsAnnoSchema(
          state.schema,
          categoryToDuplicate,
          name,
          {
            writable: true,
          }
        );
        /* if we are duplicating a non-writable annotation, it may not have an unassigned category */
        const s = schema.annotations.obsByName[name];
        if (s.categories.indexOf(unassignedCategoryLabel) === -1) {
          s.categories = s.categories.concat(unassignedCategoryLabel);
        }
        data = state.obsAnnotations.col(categoryToDuplicate).asArray();
      } else {
        /* else, all are unassined */
        const categories = [unassignedCategoryLabel];
        schema = AnnotationsHelpers.addObsAnnoSchema(state.schema, name, {
          name,
          categories,
          type: "categorical",
          writable: true,
        });
        data = new Array(state.nObs).fill(unassignedCategoryLabel);
      }

      const obsAnnotations = state.obsAnnotations.withCol(name, data);
      return { ...state, obsAnnotations, schema };
    }

    case "annotation: category edited": {
      /* change the name of an obs annotation category */
      const name = action.metadataField;
      const newName = action.newCategoryText;
      if (!AnnotationsHelpers.isUserAnnotation(state, name))
        throw new Error("unable to edit read-only annotation");
      if (typeof newName !== "string" || newName.length === 0)
        throw new Error("user annotations require string name");

      const colSchema = {
        ...state.schema.annotations.obsByName[name],
        name: newName,
      };
      const schema = AnnotationsHelpers.addObsAnnoSchema(
        AnnotationsHelpers.removeObsAnnoSchema(state.schema, name),
        newName,
        colSchema
      );
      const obsAnnotations = state.obsAnnotations.renameCol(name, newName);
      return { ...state, schema, obsAnnotations };
    }

    case "annotation: delete category": {
      /* delete annotation category from schema and obsAnnotations */
      const name = action.metadataField;
      if (!AnnotationsHelpers.isUserAnnotation(state, name))
        throw new Error("unable to delete read-only annotation");

      const schema = AnnotationsHelpers.removeObsAnnoSchema(state.schema, name);
      const obsAnnotations = state.obsAnnotations.dropCol(name);
      return { ...state, schema, obsAnnotations };
    }

    case "annotation: add new label to category": {
      const annotationName = action.metadataField;
      const newLabelName = action.newLabelText;
      if (!AnnotationsHelpers.isUserAnnotation(state, annotationName))
        throw new Error("unable to modify read-only annotation");
      if (typeof newLabelName !== "string" || newLabelName.length === 0)
        throw new Error(
          "user annotations require a non-zero length string name"
        );

      /* add the new label to the annotation schema */
      const schema = AnnotationsHelpers.addObsAnnoCategory(
        state.schema,
        annotationName,
        newLabelName
      );

      /* if so requested, label the current selection */
      const { world, crossfilter } = prevSharedState;
      const { metadataField, newLabelText } = action;
      const obsAnnotations = !action.assignSelectedCells
        ? state.obsAnnotations
        : setLabelOnCurrentSelection(
            state,
            world,
            crossfilter,
            metadataField,
            newLabelText
          );

      return { ...state, schema, obsAnnotations };
    }

    case "annotation: label edited": {
      const annotationName = action.metadataField;
      const oldLabelName = action.label;
      const newLabelName = action.editedLabel;
      if (!AnnotationsHelpers.isUserAnnotation(state, annotationName))
        throw new Error("unable to modify read-only annotation");
      if (typeof newLabelName !== "string" || newLabelName.length === 0)
        throw new Error(
          "user annotations require a non-zero length string name"
        );

      /* remove old label, add new label */
      const schema = AnnotationsHelpers.addObsAnnoCategory(
        AnnotationsHelpers.removeObsAnnoCategory(
          state.schema,
          annotationName,
          oldLabelName
        ),
        annotationName,
        newLabelName
      );

      /* change all values in obsAnnotation */
      const obsAnnotations = AnnotationsHelpers.setLabelByValue(
        state.obsAnnotations,
        annotationName,
        oldLabelName,
        newLabelName
      );

      return { ...state, schema, obsAnnotations };
    }

    case "annotation: delete label": {
      /* delete the label from the annotation, and set all cells with this value to unassigned */
      const annotationName = action.metadataField;
      const labelName = action.label;
      if (!AnnotationsHelpers.isUserAnnotation(state, annotationName))
        throw new Error("unable to modify read-only annotation");
      if (labelName === unassignedCategoryLabel)
        throw new Error("may not remove the unassigned label");

      /* remove the category from the schema */
      const schema = AnnotationsHelpers.removeObsAnnoCategory(
        state.schema,
        annotationName,
        labelName
      );

      /* set all values to unassigned in obsAnnotations */
      const obsAnnotations = AnnotationsHelpers.setLabelByValue(
        state.obsAnnotations,
        annotationName,
        labelName,
        unassignedCategoryLabel
      );

      return { ...state, schema, obsAnnotations };
    }

    case "annotation: label current cell selection": {
      const { metadataField, label } = action;
      const { world, crossfilter } = prevSharedState;
      const obsAnnotations = setLabelOnCurrentSelection(
        state,
        world,
        crossfilter,
        metadataField,
        label
      );
      return { ...state, obsAnnotations };
    }

    default: {
      return state;
    }
  }
};

function setLabelOnCurrentSelection(
  universe,
  world,
  crossfilter,
  metadataField,
  label
) {
  /*
  Set category `metadataField` to value `label` for anything currently selected.
  Used by several action type reducers.

  Returns the new obsAnnotations dataframe.
  */

  /*
  selection state is relative to world.  We need to convert it
  to a mask for Universe before applying it.
  */
  const worldMask = crossfilter.allSelectedMask();
  const mask = World.worldEqUniverse(world, universe)
    ? worldMask
    : AnnotationsHelpers.worldToUniverseMask(
        worldMask,
        world.obsAnnotations,
        universe.nObs
      );
  const obsAnnotations = AnnotationsHelpers.setLabelByMask(
    universe.obsAnnotations,
    metadataField,
    mask,
    label
  );
  return obsAnnotations;
}

export default Universe;
