/*
View creators.  These are helper functions which create new views from existing
instances of AnnoMatrix, implementing common UI functions.
*/

import { AnnoMatrixRowSubsetView, AnnoMatrixClipView } from "./views";

// eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any -- - FIXME: disabled temporarily on migrate to TS.
export function isubsetMask(annoMatrix: any, obsMask: any) {
  /*
		Subset annomatrix to contain the rows which have truish value in the mask.
    Maks length must equal annoMatrix.nObs (row count).
	*/
  return isubset(annoMatrix, _maskToList(obsMask));
}

// @ts-expect-error ts-migrate(7006) FIXME: Parameter 'annoMatrix' implicitly has an 'any' typ... Remove this comment to see the full error message
// eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types --- FIXME: disabled temporarily on migrate to TS.
export function isubset(annoMatrix, obsOffsets) {
  /*
		Subset annomatrix to contain the positions contained in the obsOffsets array

    Example:

      isubset(annoMatrix, [0, 1]) -> annoMatrix with only the first two rows
	*/
  const obsIndex = annoMatrix.rowIndex.isubset(obsOffsets);
  return new AnnoMatrixRowSubsetView(annoMatrix, obsIndex);
}

// @ts-expect-error ts-migrate(7006) FIXME: Parameter 'annoMatrix' implicitly has an 'any' typ... Remove this comment to see the full error message
// eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types --- FIXME: disabled temporarily on migrate to TS.
export function subset(annoMatrix, obsLabels) {
  /*
		subset based on labels
  */
  const obsIndex = annoMatrix.rowIndex.subset(obsLabels);
  return new AnnoMatrixRowSubsetView(annoMatrix, obsIndex);
}

// @ts-expect-error ts-migrate(7006) FIXME: Parameter 'annoMatrix' implicitly has an 'any' typ... Remove this comment to see the full error message
// eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types --- FIXME: disabled temporarily on migrate to TS.
export function subsetByIndex(annoMatrix, obsIndex) {
  /*
  subset based upon the new obs index.
  */
  return new AnnoMatrixRowSubsetView(annoMatrix, obsIndex);
}

// @ts-expect-error ts-migrate(7006) FIXME: Parameter 'annoMatrix' implicitly has an 'any' typ... Remove this comment to see the full error message
// eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types --- FIXME: disabled temporarily on migrate to TS.
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

// @ts-expect-error ts-migrate(7006) FIXME: Parameter 'mask' implicitly has an 'any' type.
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
  return list.subarray(0, elems);
}
