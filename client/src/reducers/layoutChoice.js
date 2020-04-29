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
  const idx = preferredNames.findIndex(name => layouts.indexOf(name) !== -1);
  if (idx !== -1) return preferredNames[idx];
  return layouts[0];
}

function setToDefaultLayout(world) {
  const { schema } = world;
  const available = schema.layout.obs.map(v => v.name).sort();
  const current = bestDefaultLayout(available);
  const currentDimNames = schema.layout.obsByName[current].dims;
  return { available, current, currentDimNames };
}

const LayoutChoice = (
  state = {
    available: [], // all available choices
    current: undefined, // name of the current layout, eg, 'umap'
    currentDimNames: [] // dimension name
  },
  action,
  nextSharedState
) => {
  switch (action.type) {
    case "universe exists, but loading is still in progress": {
      // set default to default
      const { universe } = nextSharedState;
      return {
        ...state,
        ...setToDefaultLayout(universe)
      };
    }

    case "set layout choice": {
      const { schema } = nextSharedState.world;
      const current = action.layoutChoice;
      const currentDimNames = schema.layout.obsByName[current].dims;
      return { ...state, current, currentDimNames };
    }

    case "reembed: add reembedding": {
      const {name} = action.schema;
      const available = Array.from(new Set(state.available).add(name));
      return {
        ...state,
        available
      };
    }

    case "reembed: clear all reembeddings": {
      const { universe } = nextSharedState;
      const { current } = state;
      const dflt = setToDefaultLayout(universe);
      if (dflt.available.includes(current)) {
        return {
          ...state,
          available: dflt.available
        };
      }
      return {
        ...state,
        ...dflt
      };
    }

    default: {
      return state;
    }
  }
};

export default LayoutChoice;
