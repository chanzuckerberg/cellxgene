import { unassignedCategoryLabel } from "../globals";
import {
  World,
  ControlsHelpers as CH,
  AnnotationsHelpers as AH
} from "../util/stateManager";
import clip from "../util/clip";
import quantile from "../util/quantile";

const WorldReducer = (
  state = null,
  action,
  nextSharedState,
  prevSharedState
) => {
  switch (action.type) {
    case "initial data load complete (universe exists)": {
      const { universe } = nextSharedState;
      const world = World.createWorldFromEntireUniverse(universe);
      return world;
    }

    case "reset World to eq Universe": {
      return prevSharedState.resetCache.world;
    }

    case "set World to current selection": {
      /* Set viewable world to be the currently selected data */
      const world = World.createWorldBySelection(
        action.universe,
        action.world,
        action.crossfilter
      );
      return world;
    }

    case "set clip quantiles": {
      const world = World.createWorldWithNewClip(
        prevSharedState.universe,
        state,
        prevSharedState.crossfilter,
        action.clipQuantiles
      );
      return world;
    }

    case "expression load success": {
      const { universe } = nextSharedState;
      const universeVarData = universe.varData;
      let unclippedVarData = state.unclipped.varData;

      // Lazy load new expression data into the unclipped varData dataframe, if
      // not already present.
      //
      Object.entries(action.expressionData).forEach(([key, val]) => {
        // If not already in world.varData, save sliced expression column
        if (!unclippedVarData.hasCol(key)) {
          // Slice if world !== universe, else just use whole column.
          // Use the obsAnnotation index as the cut key, as we keep
          // all world dataframes in sync.
          let worldValSlice = val;
          if (!World.worldEqUniverse(state, universe)) {
            worldValSlice = universeVarData
              .subset(state.obsAnnotations.rowIndex.keys(), [key], null)
              .icol(0)
              .asArray();
          }

          // Now build world's varData dataframe
          unclippedVarData = unclippedVarData.withCol(
            key,
            worldValSlice,
            state.obsAnnotations.rowIndex
          );
        }
      });

      // Prune size of varData unclipped dataframe if getting out of hand....
      //
      const { userDefinedGenes, diffexpGenes } = prevSharedState;
      const allTheGenesWeNeed = [
        ...new Set(
          userDefinedGenes,
          diffexpGenes,
          Object.keys(action.expressionData)
        )
      ];
      unclippedVarData = CH.pruneVarDataCache(
        unclippedVarData,
        allTheGenesWeNeed
      );

      // at this point, we have the unclipped data in unclippedVarData.
      // Now create clipped.
      //   - Drop columns no longer needed
      //   - Add new columns
      //
      let clippedVarData = state.varData;
      const keysToDrop = clippedVarData.colIndex
        .keys()
        .filter(k => !unclippedVarData.hasCol(k));
      const keysToAdd = unclippedVarData.colIndex
        .keys()
        .filter(k => !clippedVarData.hasCol(k));
      keysToDrop.forEach(k => {
        clippedVarData = clippedVarData.dropCol(k);
      });
      keysToAdd.forEach(k => {
        const data = unclippedVarData.col(k).asArray();
        const q = [state.clipQuantiles.min, state.clipQuantiles.max];
        const [qMinVal, qMaxVal] = quantile(q, data);
        const clippedData = clip(data, qMinVal, qMaxVal, Number.NaN);
        clippedVarData = clippedVarData.withCol(
          k,
          clippedData,
          state.obsAnnotations.rowIndex
        );
      });

      return {
        ...state,
        varData: clippedVarData,
        unclipped: {
          ...state.unclipped,
          varData: unclippedVarData
        }
      };
    }

    case "new user annotation category created": {
      const name = action.data;
      const { universe } = nextSharedState;
      const { schema } = universe;

      /*
      if world !== universe, we have to subset the newly created annotation,
      else, just use it as is.
      */
      let newAnnotation = null;
      if (!World.worldEqUniverse(state, universe)) {
        newAnnotation = universe.obsAnnotations
          .subset(state.obsAnnotations.rowIndex.keys(), [name], null)
          .icol(0)
          .asArray();
      } else {
        newAnnotation = universe.obsAnnotations.col(name).asArray();
      }
      const obsAnnotations = state.obsAnnotations.withCol(
        name,
        newAnnotation,
        state.obsAnnotations.rowIndex
      );
      return { ...state, schema, obsAnnotations };
    }

    case "category edited": {
      /* change the name of an obs annotation */
      const name = action.metadataField;
      const newName = action.editedCategoryText;
      const { schema } = nextSharedState.universe;
      const obsAnnotations = state.obsAnnotations.renameCol(name, newName);
      return { ...state, schema, obsAnnotations };
    }

    case "delete category": {
      const { schema } = nextSharedState.universe;
      const name = action.metadataField;
      const obsAnnotations = state.obsAnnotations.dropCol(name);
      return { ...state, schema, obsAnnotations };
    }

    case "add new label to category": {
      const { schema } = nextSharedState.universe;
      return { ...state, schema };
    }

    case "label edited": {
      const { schema } = nextSharedState.universe;
      const annotationName = action.metadataField;
      const oldLabelName = action.label;
      const newLabelName = action.editedLabel;

      /* set all values to to new label */
      const obsAnnotations = AH.setLabelByValue(
        state.obsAnnotations,
        annotationName,
        oldLabelName,
        newLabelName
      );
      return { ...state, schema, obsAnnotations };
    }

    case "delete label": {
      const { schema } = nextSharedState.universe;
      const annotationName = action.metadataField;
      const labelName = action.label;

      /* set all values to unassigned in obsAnnotations */
      const obsAnnotations = AH.setLabelByValue(
        state.obsAnnotations,
        annotationName,
        labelName,
        unassignedCategoryLabel
      );

      return { ...state, schema, obsAnnotations };
    }

    case "label current cell selection": {
      const { metadataField, label } = action;
      const { crossfilter } = prevSharedState;
      const mask = crossfilter.allSelectedMask();
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

export default WorldReducer;
