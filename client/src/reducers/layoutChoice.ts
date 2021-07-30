/*
we have a UI heuristic to pick the default layout, based on assumptions
about commonly used names.  Preferentially, pick in the following order:

  1. "umap"
  2. "tsne"
  3. "pca"
  4. give up, use the first available
*/
// eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
function bestDefaultLayout(layouts: any) {
  const preferredNames = ["umap", "tsne", "pca"];
  const idx = preferredNames.findIndex((name) => layouts.indexOf(name) !== -1);
  if (idx !== -1) return preferredNames[idx];
  return layouts[0];
}

// eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
function setToDefaultLayout(schema: any) {
  // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
  const available = schema.layout.obs.map((v: any) => v.name).sort();
  const current = bestDefaultLayout(available);
  const currentDimNames = schema.layout.obsByName[current].dims;
  return { available, current, currentDimNames };
}

// eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types --- FIXME: disabled temporarily on migrate to TS.
const LayoutChoice = (
  state = {
    available: [], // all available choices
    current: undefined, // name of the current layout, eg, 'umap'
    currentDimNames: [], // dimension name
  },
  // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any -- - FIXME: disabled temporarily on migrate to TS.
  action: any,
  // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any -- - FIXME: disabled temporarily on migrate to TS.
  nextSharedState: any
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
