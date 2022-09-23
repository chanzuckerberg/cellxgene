/*
View creators.  These are helper functions which create new views from existing
instances of AnnoMatrix, implementing common UI functions.
*/

import { AnnoMatrixRowSubsetView, AnnoMatrixClipView } from "./views";

export function isubsetMask(annoMatrix, obsMask) {
  /*
		Subset annomatrix to contain the rows which have truish value in the mask.
    Maks length must equal annoMatrix.nObs (row count).
	*/
  return isubset(annoMatrix, _maskToList(obsMask));
}

export function isubset(annoMatrix, obsOffsets) {
  /*
		Subset annomatrix to contain the positions contained in the obsOffsets array

    Example:

      isubset(annoMatrix, [0, 1]) -> annoMatrix with only the first two rows
	*/
  const obsIndex = annoMatrix.rowIndex.isubset(obsOffsets);
  return new AnnoMatrixRowSubsetView(annoMatrix, obsIndex);
}

export function subset(annoMatrix, obsLabels) {
  /*
		subset based on labels
  */
  const obsIndex = annoMatrix.rowIndex.subset(obsLabels);
  return new AnnoMatrixRowSubsetView(annoMatrix, obsIndex);
}

export function subsetByIndex(annoMatrix, obsIndex) {
  /*
  subset based upon the new obs index.
  */
  return new AnnoMatrixRowSubsetView(annoMatrix, obsIndex);
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
Private utility functions below
*/

function _maskToList(mask) {
  /* convert masks to lists - method wastes space, but is fast */
  if (!mask) {
    return null;
  }
  const [...m] = mask;
  const list = new Int32Array(m.length);
  let elems = 0;
  for (let i = 0, l = m.length; i < l; i += 1) {
    if (m[i]) {
      list[elems] = i;
      elems += 1;
    }
  }
  return list.subarray(0, elems);
}
