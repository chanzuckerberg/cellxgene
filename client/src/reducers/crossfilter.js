import _ from "lodash";

import Crossfilter from "../util/typedCrossfilter";
import { World, ControlsHelpers } from "../util/stateManager";
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
      const { world } = nextSharedState;
      const crossfilter = World.createObsDimensions(
        new Crossfilter(world.obsAnnotations),
        world
      );
      return crossfilter;
    }

    case "reset World to eq Universe": {
      const { userDefinedGenes, diffexpGenes } = prevSharedState.controls;
      const { world } = nextSharedState;
      const crossfilter = ControlsHelpers.createGeneDimensions(
        userDefinedGenes,
        diffexpGenes,
        world,
        prevSharedState.resetCache.crossfilter
      );
      return crossfilter;
    }

    case "set World to current selection": {
      const { userDefinedGenes, diffexpGenes } = prevSharedState.controls;
      const { world } = nextSharedState;
      let crossfilter = new Crossfilter(world.obsAnnotations);
      crossfilter = World.createObsDimensions(crossfilter, world);
      crossfilter = ControlsHelpers.createGeneDimensions(
        userDefinedGenes,
        diffexpGenes,
        world,
        crossfilter
      );
      return crossfilter;
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
      const genes = _.map(action.data, d =>
        world.varAnnotations.at(d[0], "name")
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
      const crossfilter = _.reduce(
        action.diffExp,
        (xfltr, values) => {
          const name = world.varAnnotations.at(values[0], "name");
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
      return state.select(obsAnnoDimensionName(action.metadataField), {
        mode: "exact",
        values: ControlsHelpers.selectedValuesForCategory(cat)
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
