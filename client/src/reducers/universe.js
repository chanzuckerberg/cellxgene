import { unassignedCategoryLabel } from "../globals";
import { ControlsHelpers, AnnotationsHelpers } from "../util/stateManager";

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
      varData = ControlsHelpers.pruneVarDataCache(varData, allTheGenesWeNeed);
      return { ...state, varData };
    }

    case "new user annotation category created": {
      /* create a new annotation category, with all values set to 'unassigned' */
      const name = action.data;
      /* ensure the name isn't already in use! */
      if (state.obsAnnotations.hasCol(name)) {
        console.error("name collision on annotation category create");
        return state;
      }
      const categories = [unassignedCategoryLabel];
      const schema = AnnotationsHelpers.addObsAnnoSchema(state.schema, name, {
        name,
        categories,
        type: "categorical",
        isUserAnnotation: true
      });
      const data = new Array(state.nObs).fill(unassignedCategoryLabel);
      const obsAnnotations = state.obsAnnotations.withCol(name, data);
      return { ...state, obsAnnotations, schema };
    }

    case "duplicate annotation category": {
      /* 
      create a new annotation category, with all values duplicated from an 
      existing annotation 
      */
      // TODO/XXX incomplete action - we need to know BOTH the category
      // being duplicated, and the new category name.

      return state;
    }

    case "category edited": {
      /* change the name of an obs annotation category */
      const name = action.metadataField;
      const newName = action.editedCategoryText;
      if (!AnnotationsHelpers.isUserAnnotation(state, name)) {
        console.error("deleting read-only annotation");
        return state;
      }

      const colSchema = {
        ...state.schema.annotations.obsByName[name],
        name: newName
      };
      const schema = AnnotationsHelpers.addObsAnnoSchema(
        AnnotationsHelpers.removeObsAnnoSchema(state.schema, name),
        newName,
        colSchema
      );
      const obsAnnotations = state.obsAnnotations.renameCol(name, newName);
      return { ...state, schema, obsAnnotations };
    }

    case "delete category": {
      /* delete annotation category from schema and obsAnnotations */
      const name = action.metadataField;
      if (!AnnotationsHelpers.isUserAnnotation(state, name)) {
        console.error("deleting read-only annotation");
        return state;
      }
      const schema = AnnotationsHelpers.removeObsAnnoSchema(state.schema, name);
      const obsAnnotations = state.obsAnnotations.dropCol(name);
      return { ...state, schema, obsAnnotations };
    }

    default: {
      return state;
    }
  }
};

export default Universe;
