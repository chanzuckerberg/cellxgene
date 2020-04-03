import { ColorHelpers } from "../util/stateManager";

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
    case "universe exists, but loading is still in progress": {
      /* initialize everything with default colors, no mode, no color-by accessor */
      const { world } = nextSharedState;
      const colorMode = null;
      const colorAccessor = null;
      const { rgb, scale } = ColorHelpers.createColors(world);
      return {
        ...state,
        colorAccessor,
        colorMode,
        rgb,
        scale
      };
    }

    case "universe: user color load success": {
      const { userColors } = action;
      return {
        ...state,
        userColors
      };
    }

    case "reset World to eq Universe": {
      /* need to rebuild colors as world may have changed, but don't switch modes */
      const { world } = nextSharedState;
      const { colorMode, colorAccessor } = state;
      const { rgb, scale } = ColorHelpers.createColors(
        world,
        colorMode,
        colorAccessor
      );
      return {
        ...state,
        rgb,
        scale
      };
    }

    case "set clip quantiles":
    case "set World to current selection": {
      const { world: prevWorld, controls: prevControls } = prevSharedState;
      const resetColorState = ColorHelpers.checkIfColorByDiffexpAndResetColors(
        prevControls,
        state,
        prevWorld
      );
      if (resetColorState) {
        return resetColorState;
      }

      const { colorMode, colorAccessor } = state;
      const { world } = nextSharedState;
      const { rgb, scale } = ColorHelpers.createColors(
        world,
        colorMode,
        colorAccessor
      );
      return {
        ...state,
        rgb,
        scale
      };
    }

    case "annotation: delete category": {
      const { colorAccessor } = state;
      if (action.metadataField !== colorAccessor) {
        return state;
      }
      /* else reset */
      return {
        ...state,
        ...ColorHelpers.resetColors(prevSharedState.world)
      };
    }

    case "reset colorscale": {
      return {
        ...state,
        ...ColorHelpers.resetColors(prevSharedState.world)
      };
    }

    case "color by categorical metadata":
    case "color by continuous metadata": {
      const { world } = prevSharedState;

      /* toggle between this mode and reset */
      const resetCurrent =
        action.type === state.colorMode &&
        action.colorAccessor === state.colorAccessor;
      const colorMode = !resetCurrent ? action.type : null;
      const colorAccessor = !resetCurrent ? action.colorAccessor : null;
      const { colors } = action;

      const { rgb, scale } = ColorHelpers.createColors(
        world,
        colorMode,
        colorAccessor,
        colors
      );
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
      const resetCurrent =
        action.type === state.colorMode && action.gene === state.colorAccessor;
      const colorMode = !resetCurrent ? action.type : null;
      const colorAccessor = !resetCurrent ? action.gene : null;

      const { rgb, scale } = ColorHelpers.createColors(
        world,
        colorMode,
        colorAccessor
      );
      return {
        ...state,
        colorMode,
        colorAccessor,
        rgb,
        scale
      };
    }

    case "annotation: add new label to category":
    case "annotation: label current cell selection":
    case "annotation: delete label": {
      const { world } = nextSharedState;
      const { colorMode, colorAccessor } = state;
      const { metadataField, colors } = action;
      if (
        colorMode !== "color by categorical metadata" ||
        colorAccessor !== metadataField
      )
        return state;

      /* else, we need to rebuild colors as labels have changed! */
      const { rgb, scale } = ColorHelpers.createColors(
        world,
        colorMode,
        colorAccessor,
        colors
      );
      return { ...state, rgb, scale };
    }

    case "clear differential expression": {
      const { world: prevWorld, controls: prevControls } = prevSharedState;
      const resetColorState = ColorHelpers.checkIfColorByDiffexpAndResetColors(
        prevControls,
        state,
        prevWorld
      );
      if (resetColorState) {
        return resetColorState;
      }
      return state;
    }

    default: {
      return state;
    }
  }
};

export default ColorsReducer;
