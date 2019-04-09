import { createColors } from "../util/stateManager";

const ColorsReducer = (
  state = {
    colorMode: null,
    colorAccessor: null,
    rgb: null,
    scale: null
  },
  action,
  nextSharedState,
  prevSharedState
) => {
  switch (action.type) {
    case "initial data load complete (universe exists)":
    case "reset World to eq Universe": {
      const { world } = nextSharedState;
      const colorMode = null;
      const colorAccessor = null;
      const { rgb, scale } = createColors(world, colorMode);
      return {
        ...state,
        colorAccessor,
        colorMode,
        rgb,
        scale
      };
    }

    case "set World to current selection": {
      const { colorMode, colorAccessor } = state;
      const { world } = nextSharedState;
      const { rgb, scale } = createColors(world, colorMode, colorAccessor);
      return {
        ...state,
        rgb,
        scale
      };
    }

    case "reset colorscale": {
      const { world } = prevSharedState;
      const { rgb, scale } = createColors(world);
      return {
        ...state,
        colorMode: null,
        colorAccessor: null,
        rgb,
        scale
      };
    }

    case "color by categorical metadata":
    case "color by continuous metadata": {
      const { world } = prevSharedState;

      /* toggle between this mode and reset */
      const colorMode = action.type !== state.colorMode ? action.type : null;
      const colorAccessor =
        action.colorAccessor !== state.colorAccessor && colorMode !== null
          ? action.colorAccessor
          : null;

      const { rgb, scale } = createColors(world, colorMode, colorAccessor);
      return {
        ...state,
        colorMode,
        colorAccessor,
        rgb,
        scale
      };
    }

    case "color by expression": {
      const { world } = prevSharedState;

      /* toggle between this mode and reset */
      const colorMode = action.type !== state.colorMode ? action.type : null;
      const colorAccessor =
        action.gene !== state.colorAccessor && coorMode !== null
          ? action.gene
          : null;

      const { rgb, scale } = createColors(world, colorMode, colorAccessor);
      return {
        ...state,
        colorMode,
        colorAccessor,
        rgb,
        scale
      };
    }

    default: {
      return state;
    }
  }
};

export default ColorsReducer;
