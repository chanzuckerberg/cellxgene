import { unassignedCategoryLabel } from "../globals";
import {
  World,
  ControlsHelpers,
  AnnotationsHelpers,
} from "../util/stateManager";
import {
  addObsLayout,
  removeObsLayout,
} from "../util/stateManager/schemaHelpers";
import clip from "../util/clip";
import quantile from "../util/quantile";

/*
important note: much of this code assumes that wriable (user) annotations
will NOT contain scalar data (ie, will only contain categorical labelled
data), and therefore will never need to be clipped.  Put another way, it
assumes that for these annotations, the clipped & unclipped data is equal.

If we ever start allowing user editable scalar data, this assumption will
need to be revisited.
*/

const WorldReducer = (
  state = null,
  action,
  nextSharedState,
  prevSharedState
) => {
  switch (action.type) {
    case "universe exists, but loading is still in progress":
    case "reset World to eq Universe": {
      const { universe } = nextSharedState;
      const world = World.createWorldFromEntireUniverse(universe);
      return world;
    }

    case "universe: column load success": {
      /* incremental initial data load - always assumes world == universe */
      const { universe } = nextSharedState;
      const { dim } = action;

      // we don't clip anything except for varData and obsAnnotations
      let { unclipped } = state;
      if (dim === "varData" || dim === "obsAnnotations") {
        unclipped = {
          ...unclipped,
          [dim]: universe[dim].clone(),
        };
      }
      return {
        ...state,
        schema: universe.schema,
        [dim]: universe[dim].clone(),
        unclipped,
      };
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
              .subset(state.obsAnnotations.rowIndex.labels(), [key], null)
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
          [userDefinedGenes, diffexpGenes, Object.keys(action.expressionData)]
            .filter((ele) => ele)
            .flat()
        ),
      ];
      unclippedVarData = ControlsHelpers.pruneVarDataCache(
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
        .labels()
        .filter((k) => !unclippedVarData.hasCol(k));
      const keysToAdd = unclippedVarData.colIndex
        .labels()
        .filter((k) => !clippedVarData.hasCol(k));
      keysToDrop.forEach((k) => {
        clippedVarData = clippedVarData.dropCol(k);
      });
      keysToAdd.forEach((k) => {
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
          varData: unclippedVarData,
        },
      };
    }

    case "annotation: create category": {
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
          .subset(state.obsAnnotations.rowIndex.labels(), [name], null)
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
      const unclipped = {
        ...state.unclipped,
        obsAnnotations: state.unclipped.obsAnnotations.withCol(
          name,
          newAnnotation,
          state.unclipped.obsAnnotations.rowIndex
        ),
      };
      return { ...state, schema, obsAnnotations, unclipped };
    }

    case "annotation: category edited": {
      /* change the name of an obs annotation */
      const name = action.metadataField;
      const newName = action.newCategoryText;
      const { schema } = nextSharedState.universe;
      const obsAnnotations = state.obsAnnotations.renameCol(name, newName);
      const unclipped = {
        ...state.unclipped,
        obsAnnotations: state.unclipped.obsAnnotations.renameCol(name, newName),
      };
      return { ...state, schema, obsAnnotations, unclipped };
    }

    case "annotation: delete category": {
      /* remove a category from obs annotation */
      const { schema } = nextSharedState.universe;
      const name = action.metadataField;
      const obsAnnotations = state.obsAnnotations.dropCol(name);
      const unclipped = {
        ...state.unclipped,
        obsAnnotations: state.unclipped.obsAnnotations.dropCol(name),
      };
      return { ...state, schema, obsAnnotations, unclipped };
    }

    case "annotation: add new label to category": {
      /* add a new label to the schema - schema updated by universe reducer, we just need to note it */
      const { schema } = nextSharedState.universe;
      const { metadataField, newLabelText } = action;
      const { crossfilter } = prevSharedState;

      if (action.assignSelectedCells) {
        return {
          ...state,
          schema,
          ...setLabelOnCurrentSelection(
            state,
            crossfilter,
            metadataField,
            newLabelText
          ),
        };
      }
      return { ...state, schema };
    }

    case "annotation: label edited": {
      const { schema } = nextSharedState.universe;
      const { metadataField } = action;
      const oldLabelName = action.label;
      const newLabelName = action.editedLabel;

      /* set all values to to new label */
      const unclipped = {
        ...state.unclipped,
        obsAnnotations: AnnotationsHelpers.setLabelByValue(
          state.unclipped.obsAnnotations,
          metadataField,
          oldLabelName,
          newLabelName
        ),
      };
      const obsAnnotations = state.obsAnnotations.replaceColData(
        metadataField,
        unclipped.obsAnnotations.col(metadataField).asArray()
      );
      return { ...state, schema, obsAnnotations, unclipped };
    }

    case "annotation: delete label": {
      const { schema } = nextSharedState.universe;
      const { label, metadataField } = action;

      /* set all values to unassigned in obsAnnotations */
      const unclipped = {
        ...state.unclipped,
        obsAnnotations: AnnotationsHelpers.setLabelByValue(
          state.unclipped.obsAnnotations,
          metadataField,
          label,
          unassignedCategoryLabel
        ),
      };
      const obsAnnotations = state.obsAnnotations.replaceColData(
        metadataField,
        unclipped.obsAnnotations.col(metadataField).asArray()
      );
      return { ...state, schema, obsAnnotations, unclipped };
    }

    case "annotation: label current cell selection": {
      const { metadataField, label } = action;
      const { crossfilter } = prevSharedState;
      return {
        ...state,
        ...setLabelOnCurrentSelection(state, crossfilter, metadataField, label),
      };
    }

    case "reembed: add reembedding": {
      // new embedding loaded, which *only* affects world's layout.
      // It may be new, or it may replace a previous re-embedding.
      const { obsLayout: origObsLayout, schema: origSchema } = state;
      const { embedding, schema: embeddingSchema } = action;

      const { dims } = embeddingSchema;
      let obsLayout = origObsLayout;
      let schema = origSchema;

      // alias the names the server sent us, in case they were not the same as the schema
      const embedingLabels = embedding.colIndex.labels();
      const labels = {
        [embedingLabels[0]]: dims[0],
        [embedingLabels[1]]: dims[1],
      };
      obsLayout = obsLayout.withColsFrom(embedding, labels);
      schema = addObsLayout(schema, embeddingSchema);
      return {
        ...state,
        obsLayout,
        schema,
      };
    }

    case "reembed: clear all reembeddings": {
      // reembedding was cleared -- remove from layout
      const { obsLayout: origObsLayout, schema: origSchema } = state;
      const { reembedding } = prevSharedState;

      let schema = origSchema;
      let obsLayout = origObsLayout;

      reembedding.reembeddings.forEach((emb, name) => {
        const { dims } = emb.schema;
        obsLayout = obsLayout.dropCol(dims[0]);
        obsLayout = obsLayout.dropCol(dims[1]);
        schema = removeObsLayout(schema, name);
      });
      return {
        ...state,
        obsLayout,
        schema,
      };
    }

    default: {
      return state;
    }
  }
};

function setLabelOnCurrentSelection(world, crossfilter, metadataField, label) {
  /*
  Set category `metadataField` to value `label` for anything currently selected.
  Used by several action type reducers.
  */
  const mask = crossfilter.allSelectedMask();
  const unclipped = {
    ...world.unclipped,
    obsAnnotations: AnnotationsHelpers.setLabelByMask(
      world.unclipped.obsAnnotations,
      metadataField,
      mask,
      label
    ),
  };
  const obsAnnotations = world.obsAnnotations.replaceColData(
    metadataField,
    unclipped.obsAnnotations.col(metadataField).asArray()
  );

  return { obsAnnotations, unclipped };
}

export default WorldReducer;
