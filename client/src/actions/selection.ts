/*
Action creators for selection 
*/
export const selectContinuousMetadataAction = (
  type: any,
  query: any,
  range: any,
  oldProps = {}
) => async (dispatch: any, getState: any) => {
  const { obsCrossfilter: prevObsCrossfilter } = getState();

  const selection = range
    ? {
        mode: "range",
        lo: range[0],
        hi: range[1],
        inclusive: true, // [lo, hi] incluisve selection
      }
    : { mode: "all" };

  const obsCrossfilter = await prevObsCrossfilter.select(...query, selection);

  dispatch({
    type,
    obsCrossfilter,
    range,
    ...oldProps,
  });
};

export const selectCategoricalMetadataAction = (
  type: any, // action type
  metadataField: any, // annotation category name
  labels: any,
  label: any, // the label being selected/deselected
  isSelected: any, // bool
  oldProps = {}
) => async (dispatch: any, getState: any) => {
  const {
    obsCrossfilter: prevObsCrossfilter,
    categoricalSelection,
  } = getState();

  const labelSelectionState = new Map(categoricalSelection[metadataField]);
  labels.forEach(
    (l: any) => labelSelectionState.has(l) || labelSelectionState.set(l, true)
  );
  labelSelectionState.set(label, isSelected);

  const values = Array.from(labelSelectionState.keys()).filter((k) =>
    labelSelectionState.get(k)
  );
  const selection = {
    mode: "exact",
    values,
  };
  const obsCrossfilter = await prevObsCrossfilter.select(
    "obs",
    metadataField,
    selection
  );

  dispatch({
    type,
    obsCrossfilter,
    metadataField,
    labelSelectionState,
    ...oldProps,
  });
};

export const selectCategoricalAllMetadataAction = (
  type: any, // action type
  metadataField: any, // annotation category name
  labels: any,
  isSelected: any, // bool, select all or none
  oldProps = {}
) => async (dispatch: any, getState: any) => {
  const {
    obsCrossfilter: prevObsCrossfilter,
    categoricalSelection,
  } = getState();

  const labelSelectionState = new Map(categoricalSelection[metadataField]);
  labels.forEach((label: any) => labelSelectionState.set(label, isSelected));

  const selection = { mode: isSelected ? "all" : "none" };
  const obsCrossfilter = await prevObsCrossfilter.select(
    "obs",
    metadataField,
    selection
  );

  dispatch({
    type,
    obsCrossfilter,
    metadataField,
    labelSelectionState,
    ...oldProps,
  });
};

/**
 ** Graph selection-related actions
 **/

export const graphBrushStartAction = () =>
  /* no change to crossfilter until a change fires */
  ({ type: "graph brush start" });

const _graphBrushWithinRectAction = (
  type: any,
  embName: any,
  brushCoords: any
) => async (dispatch: any, getState: any) => {
  const { obsCrossfilter: prevObsCrossfilter } = getState();

  const selection = { mode: "within-rect", ...brushCoords };
  const obsCrossfilter = await prevObsCrossfilter.select(
    "emb",
    embName,
    selection
  );

  dispatch({
    type,
    obsCrossfilter,
    brushCoords,
  });
};

const _graphAllAction = (type: any, embName: any) => async (
  dispatch: any,
  getState: any
) => {
  const { obsCrossfilter: prevObsCrossfilter } = getState();

  const obsCrossfilter = await prevObsCrossfilter.select("emb", embName, {
    mode: "all",
  });

  dispatch({
    type,
    obsCrossfilter,
  });
};

export const graphBrushChangeAction = (embName: any, brushCoords: any) =>
  _graphBrushWithinRectAction("graph brush change", embName, brushCoords);

export const graphBrushEndAction = (embName: any, brushCoords: any) =>
  _graphBrushWithinRectAction("graph brush end", embName, brushCoords);

export const graphBrushCancelAction = (embName: any) =>
  _graphAllAction("graph brush cancel", embName);
export const graphBrushDeselectAction = (embName: any) =>
  _graphAllAction("graph brush deselect", embName);

export const graphLassoStartAction = () =>
  /* no change to crossfilter until a change fires */
  ({ type: "graph lasso start" });

export const graphLassoCancelAction = (embName: any) =>
  _graphAllAction("graph lasso cancel", embName);

export const graphLassoDeselectAction = (embName: any) =>
  _graphAllAction("graph lasso cancel", embName);

export const graphLassoEndAction = (embName: any, polygon: any) => async (
  dispatch: any,
  getState: any
) => {
  const { obsCrossfilter: prevObsCrossfilter } = getState();

  const selection = {
    mode: "within-polygon",
    polygon,
  };
  const obsCrossfilter = await prevObsCrossfilter.select(
    "emb",
    embName,
    selection
  );

  dispatch({
    type: "graph lasso end",
    obsCrossfilter,
    polygon,
  });
};

/*
Differential expression set selection
*/
export const setCellSetFromSelection = (cellSetId: any) => (
  dispatch: any,
  getState: any
) => {
  const { obsCrossfilter } = getState();
  const selected = obsCrossfilter.allSelectedLabels();

  dispatch({
    type: `store current cell selection as differential set ${cellSetId}`,
    data: selected.length > 0 ? selected : null,
  });
};
