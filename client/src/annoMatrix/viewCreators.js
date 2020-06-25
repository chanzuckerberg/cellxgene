/*
View creators
*/

import { AnnoMatrixRowSubsetView, AnnoMatrixClipView } from "./views";

export function isubsetMask(annoMatrix, rowMask) {
  /*
		Subset on row based upon mask
		*/
  return isubset(annoMatrix, _maskToList(rowMask));
}

export function isubset(annoMatrix, rowOffsets) {
  /*
		subset based on offset
		*/
  const rowIndex = annoMatrix.rowIndex.isubset(rowOffsets);
  return new AnnoMatrixRowSubsetView(annoMatrix, rowIndex);
}

export function subset(annoMatrix, rowLabels) {
  /*
		subset based on labels
		*/
  const rowIndex = annoMatrix.rowIndex.subset(rowLabels);
  return new AnnoMatrixRowSubsetView(annoMatrix, rowIndex);
}

export function clip(annoMatrix, qmin, qmax) {
  /*
		Create a view that clips all continuous data to the [min, max] range.
		The matrix shape does not change, but the continuous values outside the
		specified range will become a NaN.
		*/
  return new AnnoMatrixClipView(annoMatrix, qmin, qmax);
}

/*
Utility functions below
*/

function _maskToList(mask) {
  /* convert masks to lists - method wastes space, but is fast */
  if (!mask) {
    return null;
  }
  const list = new Int32Array(mask.length);
  let elems = 0;
  for (let i = 0, l = mask.length; i < l; i += 1) {
    if (mask[i]) {
      list[elems] = i;
      elems += 1;
    }
  }
  return new Int32Array(list.buffer, 0, elems);
}
