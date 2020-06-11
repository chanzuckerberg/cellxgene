import _ from "lodash";

import Crossfilter from "../util/typedCrossfilter";
import {
  World,
  ControlsHelpers,
  AnnotationsHelpers,
} from "../util/stateManager";
import {
  layoutDimensionName,
  obsAnnoDimensionName,
  userDefinedDimensionName,
  diffexpDimensionName,
  makeContinuousDimensionName,
} from "../util/nameCreators";

const XYDimName = layoutDimensionName("XY");

const CrossfilterReducerBase = (
  state = null,
  action,
  nextSharedState,
  prevSharedState
) => {
  switch (action.type) {
    case "universe: column load success": {
      const { world, layoutChoice } = nextSharedState;
      const { obsAnnotations, obsLayout } = world;

      // ignore var dimension loads as these are not currently selectable
      if (action.dim === "varAnnotations") return state;

      /* 
      during bootstrap loading, we don't know if obsLayout or obsAnnotations
      will load first. Take whichever arrives and is not empty (so that our
      crossfilter has the right dimensionality).
      */
      let crossfilter =
        state ??
        new Crossfilter(obsAnnotations.isEmpty() ? obsLayout : obsAnnotations);

      // add layout dimension, if not already present
      if (
        obsLayout.hasCol(layoutChoice.currentDimNames[0]) &&
        !crossfilter.hasDimension(XYDimName)
      ) {
        crossfilter = crossfilter.addDimension(
          XYDimName,
          "spatial",
          obsLayout.col(layoutChoice.currentDimNames[0]).asArray(),
          obsLayout.col(layoutChoice.currentDimNames[1]).asArray()
        );
      }

      // add any missing obsAnnotations
      crossfilter = World.addObsDimensions(crossfilter, world);
      return crossfilter;
    }

    case "reset World to eq Universe": {
      const { userDefinedGenes, diffexpGenes } = nextSharedState.controls;
      const { world } = nextSharedState;
      let { crossfilter } = prevSharedState.resetCache;
      crossfilter = ControlsHelpers.createGeneDimensions(
        userDefinedGenes,
        diffexpGenes,
        world,
        crossfilter
      );
      crossfilter = AnnotationsHelpers.createWritableAnnotationDimensions(
        world,
        crossfilter
      );
      return crossfilter;
    }

    case "set clip quantiles":
    case "set World to current selection": {
      const { userDefinedGenes, diffexpGenes } = nextSharedState.controls;
      const { world, layoutChoice } = nextSharedState;
      let crossfilter = new Crossfilter(world.obsAnnotations);
      crossfilter = World.createObsDimensions(
        crossfilter,
        world,
        layoutChoice.currentDimNames
      );
      crossfilter = ControlsHelpers.createGeneDimensions(
        userDefinedGenes,
        diffexpGenes,
        world,
        crossfilter
      );
      return crossfilter;
    }

    case "set layout choice": {
      /*
      when switching layouts:
      - delete the existing XY index
      - add the new XY index (which implicitly selects all on it)
      */
      const { world, layoutChoice } = nextSharedState;
      return state
        .delDimension(layoutDimensionName("XY"))
        .addDimension(
          layoutDimensionName("XY"),
          "spatial",
          world.obsLayout.col(layoutChoice.currentDimNames[0]).asArray(),
          world.obsLayout.col(layoutChoice.currentDimNames[1]).asArray()
        );
    }

    case "request user defined gene success": {
      const { world } = nextSharedState;
      const gene = action.data.genes[0];
      return state.addDimension(
        userDefinedDimensionName(gene),
        "scalar",
        world.varData.col(gene).asArray(),
        Float32Array
      );
    }

    case "request differential expression success": {
      const { world } = nextSharedState;
      const varIndexName = world.schema.annotations.var.index;
      const genes = _.map(action.data, (d) =>
        world.varAnnotations.at(d[0], varIndexName)
      );
      const crossfilter = _.reduce(
        genes,
        (xfltr, gene) =>
          xfltr.addDimension(
            diffexpDimensionName(gene),
            "scalar",
            world.varData.col(gene).asArray(),
            Float32Array
          ),
        state
      );
      return crossfilter;
    }

    case "clear differential expression": {
      const { world } = nextSharedState;
      const varIndexName = world.schema.annotations.var.index;
      const crossfilter = _.reduce(
        action.diffExp,
        (xfltr, values) => {
          const name = world.varAnnotations.at(values[0], varIndexName);
          return xfltr.delDimension(diffexpDimensionName(name));
        },
        state
      );
      return crossfilter;
    }

    case "clear user defined gene": {
      return state.delDimension(userDefinedDimensionName(action.data));
    }

    case "clear all user defined genes": {
      const { userDefinedGenes } = nextSharedState.controls;
      const crossfilter = _.reduce(
        userDefinedGenes,
        (xfltr, gene) => xfltr.delDimension(userDefinedDimensionName(gene)),
        state
      );
      return crossfilter;
    }

    case "annotation: create category": {
      const name = action.data;
      const { world } = nextSharedState;
      const colData = world.obsAnnotations.col(name).asArray();
      return state.addDimension(obsAnnoDimensionName(name), "enum", colData);
    }

    case "annotation: category edited": {
      const name = action.metadataField;
      const newName = action.newCategoryText;
      return state.renameDimension(
        obsAnnoDimensionName(name),
        obsAnnoDimensionName(newName)
      );
    }

    case "annotation: delete category": {
      return state.delDimension(obsAnnoDimensionName(action.metadataField));
    }

    case "annotation: add new label to category":
    case "annotation: label current cell selection":
    case "annotation: label edited":
    case "annotation: delete label": {
      if (
        action.type === "annotation: add new label to category" &&
        !action.assignSelectedCells
      )
        return state;

      /* we need to reindex the dimension.  For now, just drop it and add another */
      const name = action.metadataField;
      const dimName = obsAnnoDimensionName(name);
      const { world } = nextSharedState;
      const colData = world.obsAnnotations.col(name).asArray();
      return state.delDimension(dimName).addDimension(dimName, "enum", colData);
    }

    case "graph brush end":
    case "graph brush change": {
      const [minX, maxY] = action.brushCoords.northwest;
      const [maxX, minY] = action.brushCoords.southeast;
      return state.select(XYDimName, {
        mode: "within-rect",
        minX,
        minY,
        maxX,
        maxY,
      });
    }

    case "graph lasso end": {
      const { polygon } = action;
      return state.select(XYDimName, {
        mode: "within-polygon",
        polygon,
      });
    }

    case "graph lasso cancel":
    case "graph brush cancel":
    case "graph lasso deselect":
    case "graph brush deselect": {
      return state.select(XYDimName, { mode: "all" });
    }

    case "continuous metadata histogram start":
    case "continuous metadata histogram brush":
    case "continuous metadata histogram cancel":
    case "continuous metadata histogram end": {
      const name = makeContinuousDimensionName(
        action.continuousNamespace,
        action.selection
      );
      // action.selection: metadata name being selected
      // action.range: filter range, or null if deselected
      if (!action.range) {
        return state.select(name, { mode: "all" });
      }
      const [lo, hi] = action.range;
      const newState = state.select(name, {
        mode: "range",
        lo,
        hi,
        inclusive: true, // [lo, hi] incluisve selection
      });
      return newState;
    }

    case "categorical metadata filter select":
    case "categorical metadata filter deselect": {
      const { labels, metadataField } = action;
      const selected = nextSharedState.categoricalSelection[metadataField];
      const values = labels.filter((label) => selected.get(label) ?? true);
      return state.select(obsAnnoDimensionName(action.metadataField), {
        mode: "exact",
        values,
      });
    }

    case "categorical metadata filter none of these": {
      return state.select(obsAnnoDimensionName(action.metadataField), {
        mode: "none",
      });
    }

    case "categorical metadata filter all of these": {
      return state.select(obsAnnoDimensionName(action.metadataField), {
        mode: "all",
      });
    }

    default: {
      return state;
    }
  }
};

/*
  IMPORTANT: the system assumes that crossfilter.data() will point at the
  same value as world.obsAnnotations.  For actions handled in this reducer,
  make sure that this remains true.

  This wrapper performs only this function.
*/
const CrossfilterReducer = (
  state,
  action,
  nextSharedState,
  prevSharedState
) => {
  const nextState = CrossfilterReducerBase(
    state,
    action,
    nextSharedState,
    prevSharedState
  );
  /*
  update the data in the crossfilter to point at the current obsAnnotations, IF
  they are not empty.  If empty, leave it alone (can occur during boostrap loading).
  */
  const nextObsAnnotations = nextSharedState.world?.obsAnnotations;
  if (
    !nextState ||
    nextState.all() === nextObsAnnotations ||
    nextObsAnnotations.isEmpty()
  ) {
    return nextState;
  }
  return nextState.setData(nextObsAnnotations);
};

export default CrossfilterReducer;
