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
  let obsCrossfilter;
  if (type === "graph lasso deselect") {
    obsCrossfilter = await prevObsCrossfilter.select("obs", "name_0", {
      mode: "all",
    });
  } else {
    obsCrossfilter = prevObsCrossfilter;
  }
  obsCrossfilter = await obsCrossfilter.select("emb", embName, {
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
  _graphAllAction("graph lasso deselect", embName);

export const graphLassoEndAction = (embName, polygon, multiselect) => async (
  dispatch,
  getState
) => {
  const { annoMatrix, obsCrossfilter: prevObsCrossfilter } = getState();
  const isMultiselectOn = multiselect ?? false;
  let selection = {
    mode: "within-polygon",
    polygon,
  };
  let obsCrossfilter = await prevObsCrossfilter.select("obs", "name_0", {
    mode: "all",
  });
  obsCrossfilter = await obsCrossfilter.select("emb", embName, selection);

  const embSelection =
    prevObsCrossfilter.obsCrossfilter.dimensions[
      `emb/${embName}_0:${embName}_1`
    ]?.selection;
  let countEmb = 0;
  if (embSelection) {
    const { ranges } = embSelection;
    ranges.forEach((range) => {
      countEmb += range[1] - range[0];
    });
  } else {
    countEmb = annoMatrix.nObs;
  }

  const nameSelection =
    prevObsCrossfilter.obsCrossfilter.dimensions[`obs/name_0`]?.selection;
  let countName = 0;
  if (nameSelection) {
    const { ranges } = nameSelection;
    ranges.forEach((range) => {
      countName += range[1] - range[0];
    });
  } else {
    countName = annoMatrix.nObs;
  }

  if (
    isMultiselectOn &&
    (countEmb < annoMatrix.nObs || countName < annoMatrix.nObs)
  ) {
    const arr1 = new Array(annoMatrix.nObs);
    const arr2 = new Array(annoMatrix.nObs);
    obsCrossfilter.fillByIsSelected(arr1, true, false);
    prevObsCrossfilter.fillByIsSelected(arr2, true, false);
    const union = [...arr1.keys()].filter((i) => arr1[i] || arr2[i]);
    const nameDf = await annoMatrix.fetch("obs", "name_0");
    const rowNames = nameDf.__columns[0];
    const values = union.map((index) => rowNames[index]);

    selection = {
      mode: "exact",
      values,
    };
    obsCrossfilter = await obsCrossfilter.select("obs", "name_0", selection);
    obsCrossfilter = await obsCrossfilter.select("emb", embName, {
      mode: "all",
    });
  }

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
