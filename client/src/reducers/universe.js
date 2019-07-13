import { unassignedCategoryLabel } from "../globals";
import {
  World,
  ControlsHelpers as CH,
  AnnotationsHelpers as AH
} from "../util/stateManager";

const Universe = (state = null, action, nextSharedState, prevSharedState) => {
  switch (action.type) {
    case "initial data load complete (universe exists)": {
      const { universe } = action;
      return universe;
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
          userDefinedGenes,
          diffexpGenes,
          Object.keys(action.expressionData)
        )
      ];
      varData = CH.pruneVarDataCache(varData, allTheGenesWeNeed);
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
        schema = AH.dupObsAnnoSchema(state.schema, categoryToDuplicate, name, {
          writable: true
        });
        data = state.obsAnnotations.col(categoryToDuplicate).asArray();
      } else {
        /* else, all are unassined */
        const categories = [unassignedCategoryLabel];
        schema = AH.addObsAnnoSchema(state.schema, name, {
          name,
          categories,
          type: "categorical",
          writable: true
        });
        data = new Array(state.nObs).fill(unassignedCategoryLabel);
      }

      const obsAnnotations = state.obsAnnotations.withCol(name, data);
      return { ...state, obsAnnotations, schema };
    }

    case "annotation: category edited": {
      /* change the name of an obs annotation category */
      const name = action.metadataField;
      const newName = action.editedCategoryText;
      if (!AH.isUserAnnotation(state, name))
        throw new Error("unable to edit read-only annotation");
      if (typeof newName !== "string" || newName.length === 0)
        throw new Error("user annotations require string name");

      const colSchema = {
        ...state.schema.annotations.obsByName[name],
        name: newName
      };
      const schema = AH.addObsAnnoSchema(
        AH.removeObsAnnoSchema(state.schema, name),
        newName,
        colSchema
      );
      const obsAnnotations = state.obsAnnotations.renameCol(name, newName);
      return { ...state, schema, obsAnnotations };
    }

    case "annotation: delete category": {
      /* delete annotation category from schema and obsAnnotations */
      const name = action.metadataField;
      if (!AH.isUserAnnotation(state, name))
        throw new Error("unable to delete read-only annotation");

      const schema = AH.removeObsAnnoSchema(state.schema, name);
      const obsAnnotations = state.obsAnnotations.dropCol(name);
      return { ...state, schema, obsAnnotations };
    }

    case "annotation: add new label to category": {
      const annotationName = action.metadataField;
      const newLabelName = action.newLabelText;
      if (!AH.isUserAnnotation(state, annotationName))
        throw new Error("unable to modify read-only annotation");
      if (typeof newLabelName !== "string" || newLabelName.length === 0)
        throw new Error(
          "user annotations require a non-zero length string name"
        );

      /* add the new label to the annotation */
      const schema = AH.addObsAnnoCategory(
        state.schema,
        annotationName,
        newLabelName
      );
      return { ...state, schema };
    }

    case "annotation: label edited": {
      const annotationName = action.metadataField;
      const oldLabelName = action.label;
      const newLabelName = action.editedLabel;
      if (!AH.isUserAnnotation(state, annotationName))
        throw new Error("unable to modify read-only annotation");
      if (typeof newLabelName !== "string" || newLabelName.length === 0)
        throw new Error(
          "user annotations require a non-zero length string name"
        );

      /* remove old label, add new label */
      const schema = AH.addObsAnnoCategory(
        AH.removeObsAnnoCategory(state.schema, annotationName, oldLabelName),
        annotationName,
        newLabelName
      );

      /* change all values in obsAnnotation */
      const obsAnnotations = AH.setLabelByValue(
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
      if (!AH.isUserAnnotation(state, annotationName))
        throw new Error("unable to modify read-only annotation");
      if (labelName === unassignedCategoryLabel)
        throw new Error("may not remove the unassigned label");

      /* remove the category from the schema */
      const schema = AH.removeObsAnnoCategory(
        state.schema,
        annotationName,
        labelName
      );

      /* set all values to unassigned in obsAnnotations */
      const obsAnnotations = AH.setLabelByValue(
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

      /*
      selection state is relative to world.  We need to convert it
      to a mask for Universe before applying it.
      */
      const worldMask = crossfilter.allSelectedMask();
      const mask = World.worldEqUniverse(world, state)
        ? worldMask
        : AH.worldToUniverseMask(worldMask, world.obsAnnotations, state.nObs);
      const obsAnnotations = AH.setLabelByMask(
        state.obsAnnotations,
        metadataField,
        mask,
        label
      );
      return { ...state, obsAnnotations };
    }

    default: {
      return state;
    }
  }
};

export default Universe;
