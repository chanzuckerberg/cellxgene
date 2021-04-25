/*
we have a UI heuristic to pick the default layout, based on assumptions
about commonly used names.  Preferentially, pick in the following order:

  1. "umap"
  2. "tsne"
  3. "pca"
  4. give up, use the first available
*/
function bestDefaultLayout(layouts) {
  const preferredNames = ["umap", "tsne", "pca"];
  const idx = preferredNames.findIndex((name) => layouts.indexOf(name) !== -1);
  if (idx !== -1) return preferredNames[idx];
  return layouts[0];
}

function setToDefaultLayout(schema) {
  const available = schema.layout.obs.map((v) => v.name).sort();
  const current = bestDefaultLayout(available);
  const currentDimNames = schema.layout.obsByName[current].dims;
  return { available, current, currentDimNames };
}

const LayoutChoice = (
  state = {
    available: [], // all available choices
    dotplot: true /* todo(colinmegill) #632 wire up inputs */,
    current: undefined, // name of the current layout, eg, 'umap'
    currentDimNames: [], // dimension name
  },
  action,
  nextSharedState
) => {
  switch (action.type) {
    case "initial data load complete": {
      // set default to default
      const { annoMatrix } = nextSharedState;
      return {
        ...state,
        ...setToDefaultLayout(annoMatrix.schema),
      };
    }

    case "set layout choice": {
      const { schema } = nextSharedState.annoMatrix;
      const current = action.layoutChoice;
      const currentDimNames = schema.layout.obsByName[current].dims;
      return { ...state, current, currentDimNames };
    }

    default: {
      return state;
    }
  }
};

export default LayoutChoice;
