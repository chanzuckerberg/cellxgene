/*
Action creators for selection 
*/
export const selectContinuousMetadataAction = (
  type,
  query,
  range,
  oldProps = {}
) => async (dispatch, getState) => {
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
  type, // action type
  metadataField, // annotation category name
  labels,
  label, // the label being selected/deselected
  isSelected, // bool
  oldProps = {}
) => async (dispatch, getState) => {
  const {
    obsCrossfilter: prevObsCrossfilter,
    categoricalSelection,
  } = getState();

  const labelSelectionState = new Map(categoricalSelection[metadataField]);
  labels.forEach(
    (l) => labelSelectionState.has(l) || labelSelectionState.set(l, true)
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
  type, // action type
  metadataField, // annotation category name
  labels,
  isSelected, // bool, select all or none
  oldProps = {}
) => async (dispatch, getState) => {
  const {
    obsCrossfilter: prevObsCrossfilter,
    categoricalSelection,
  } = getState();

  const labelSelectionState = new Map(categoricalSelection[metadataField]);
  labels.forEach((label) => labelSelectionState.set(label, isSelected));

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

const _graphBrushWithinRectAction = (type, embName, brushCoords) => async (
  dispatch,
  getState
) => {
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

const _graphAllAction = (type, embName) => async (dispatch, getState) => {
  const { obsCrossfilter: prevObsCrossfilter } = getState();

  const obsCrossfilter = await prevObsCrossfilter.select("emb", embName, {
    mode: "all",
  });

  dispatch({
    type,
    obsCrossfilter,
  });
};

export const graphBrushChangeAction = (embName, brushCoords) =>
  _graphBrushWithinRectAction("graph brush change", embName, brushCoords);

export const graphBrushEndAction = (embName, brushCoords) =>
  _graphBrushWithinRectAction("graph brush end", embName, brushCoords);

export const graphBrushCancelAction = (embName) =>
  _graphAllAction("graph brush cancel", embName);
export const graphBrushDeselectAction = (embName) =>
  _graphAllAction("graph brush deselect", embName);

export const graphLassoStartAction = () =>
  /* no change to crossfilter until a change fires */
  ({ type: "graph lasso start" });

export const graphLassoCancelAction = (embName) =>
  _graphAllAction("graph lasso cancel", embName);

export const graphLassoDeselectAction = (embName) =>
  _graphAllAction("graph lasso cancel", embName);

export const graphLassoEndAction = (embName, polygon) => async (
  dispatch,
  getState
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

export const setCellsFromSelectionAndInverseAction = () => async (
  dispatch,
  getState
) => {
  const { annoMatrix, obsCrossfilter: prevObsCrossfilter } = getState();
  const arr = new Array(annoMatrix.nObs);
  prevObsCrossfilter.fillByIsSelected(arr, false, true);

  const unselected = [...arr.keys()].filter((i) => arr[i]);
  dispatch(setCellSetFromSelection(1));
  dispatch(setCellSetFromInputArray(2, unselected));
};

/*
Differential expression set selection
*/
export const setCellSetFromSelection = (cellSetId) => (dispatch, getState) => {
  const { obsCrossfilter } = getState();
  const selected = obsCrossfilter.allSelectedLabels();

  dispatch({
    type: `store current cell selection as differential set ${cellSetId}`,
    data: selected.length > 0 ? selected : null,
  });
};

export const setCellSetFromInputArray = (cellSetId, cells) => (dispatch) => {
  dispatch({
    type: `store current cell selection as differential set ${cellSetId}`,
    data: cells.length > 0 ? cells : null,
  });
};
