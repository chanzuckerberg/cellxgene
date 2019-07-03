import _ from "lodash";

import Crossfilter from "../util/typedCrossfilter";
import { World, ControlsHelpers as CH } from "../util/stateManager";
import {
  layoutDimensionName,
  obsAnnoDimensionName,
  userDefinedDimensionName,
  diffexpDimensionName,
  makeContinuousDimensionName
} from "../util/nameCreators";

const XYDimName = layoutDimensionName("XY");

const CrossfilterReducer = (
  state = null,
  action,
  nextSharedState,
  prevSharedState
) => {
  switch (action.type) {
    case "initial data load complete (universe exists)": {
      const { world, layoutChoice } = nextSharedState;
      const crossfilter = World.createObsDimensions(
        new Crossfilter(world.obsAnnotations),
        world,
        layoutChoice.currentDimNames
      );
      return crossfilter;
    }

    case "reset World to eq Universe": {
      const { userDefinedGenes, diffexpGenes } = prevSharedState.controls;
      const { world } = nextSharedState;
      const crossfilter = CH.createGeneDimensions(
        userDefinedGenes,
        diffexpGenes,
        world,
        prevSharedState.resetCache.crossfilter
      );
      return crossfilter;
    }

    case "set clip quantiles":
    case "set World to current selection": {
      const { userDefinedGenes, diffexpGenes } = prevSharedState.controls;
      const { world, layoutChoice } = nextSharedState;
      let crossfilter = new Crossfilter(world.obsAnnotations);
      crossfilter = World.createObsDimensions(
        crossfilter,
        world,
        layoutChoice.currentDimNames
      );
      crossfilter = CH.createGeneDimensions(
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
      const { world } = prevSharedState;
      const gene = action.data.genes[0];
      return state.addDimension(
        userDefinedDimensionName(gene),
        "scalar",
        world.varData.col(gene).asArray(),
        Float32Array
      );
    }

    case "request differential expression success": {
      const { world } = prevSharedState;
      const varIndexName = world.schema.annotations.var.index;
      const genes = _.map(action.data, d =>
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
      const { world } = prevSharedState;
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
      const { userDefinedGenes } = prevSharedState.controls;
      const crossfilter = _.reduce(
        userDefinedGenes,
        (xfltr, gene) => xfltr.delDimension(userDefinedDimensionName(gene)),
        state
      );
      return crossfilter;
    }

    case "new user annotation category created": {
      const name = action.data;
      const { world } = nextSharedState;
      const colData = world.obsAnnotations.col(name).asArray();
      return state.addDimension(obsAnnoDimensionName(name), "enum", colData);
    }

    case "duplicate annotation category": {
      // TODO: this code is right, but waiting on world/universe reducers
      // const name = action.metadataField;
      // const { world } = nextSharedState;
      // const colData = world.obsAnnotations.col(name).asArray();
      // return state.addDimension(obsAnnoDimensionName(name), "enum", colData);
      return state;
    }

    case "category edited": {
      const name = action.metadataField;
      const newName = action.editedCategoryText;
      return state.renameDimension(
        obsAnnoDimensionName(name),
        obsAnnoDimensionName(newName)
      );
    }

    case "delete category": {
      return state.delDimension(obsAnnoDimensionName(action.metadataField));
    }

    case "label current cell selection":
    case "label edited":
    case "delete label": {
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
        maxY
      });
    }

    case "graph lasso end": {
      const { polygon } = action;
      return state.select(XYDimName, {
        mode: "within-polygon",
        polygon
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
      const newState = state.select(name, { mode: "range", lo, hi });
      return newState;
    }

    case "categorical metadata filter select":
    case "categorical metadata filter deselect": {
      const { categoricalSelection } = nextSharedState;
      const cat = categoricalSelection[action.metadataField];
      const { categoryValues, categoryValueSelected } = cat;
      const values = categoryValues.filter((v, i) => categoryValueSelected[i]);
      return state.select(obsAnnoDimensionName(action.metadataField), {
        mode: "exact",
        values
      });
    }

    case "categorical metadata filter none of these": {
      return state.select(obsAnnoDimensionName(action.metadataField), {
        mode: "none"
      });
    }

    case "categorical metadata filter all of these": {
      return state.select(obsAnnoDimensionName(action.metadataField), {
        mode: "all"
      });
    }

    default: {
      return state;
    }
  }
};

export default CrossfilterReducer;
