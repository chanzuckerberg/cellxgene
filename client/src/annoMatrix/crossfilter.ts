/*
Row crossfilter proxy for an AnnoMatrix.  This wraps Crossfilter,
providing a number of services, and ensuring that the crossfilter and
AnnoMatrix stay in sync:
	- on-demand index creation as data is loaded
	- transparently mapping between queries and crossfilter index names.
  - for mutation of the matrix by user annotations, maintain synchronization
    between Crossfilter and AnnoMatrix.
*/
import Crossfilter from "../util/typedCrossfilter";
import { _getColumnSchema } from "./schema";

// @ts-expect-error ts-migrate(7006) FIXME: Parameter 'field' implicitly has an 'any' type.
function _dimensionNameFromDf(field, df) {
  const colNames = df.colIndex.labels();
  return _dimensionName(field, colNames);
}

// @ts-expect-error ts-migrate(7006) FIXME: Parameter 'field' implicitly has an 'any' type.
function _dimensionName(field, colNames) {
  if (!Array.isArray(colNames)) return `${field}/${colNames}`;
  return `${field}/${colNames.join(":")}`;
}

export default class AnnoMatrixObsCrossfilter {
  // @ts-expect-error ts-migrate(7006) FIXME: Parameter 'annoMatrix' implicitly has an 'any' typ... Remove this comment to see the full error message
  // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types --- FIXME: disabled temporarily on migrate to TS.
  constructor(annoMatrix, _obsCrossfilter = null) {
    // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
    (this as any).annoMatrix = annoMatrix;
    // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
    (this as any).obsCrossfilter =
      _obsCrossfilter || new Crossfilter(annoMatrix._cache.obs);
    // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
    (this as any).obsCrossfilter = (this as any).obsCrossfilter.setData(
      annoMatrix._cache.obs
    );
  }

  // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types --- FIXME: disabled temporarily on migrate to TS.
  size() {
    // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
    return (this as any).obsCrossfilter.size();
  }

  /**
  Managing the associated annoMatrix.  These wrappers are necessary to
  make coordinated changes to BOTH the crossfilter and annoMatrix, and
  ensure that all state stays synchronized.

  See API documentation in annoMatrix.js.
  **/
  // @ts-expect-error ts-migrate(7006) FIXME: Parameter 'colSchema' implicitly has an 'any' type... Remove this comment to see the full error message
  // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types --- FIXME: disabled temporarily on migrate to TS.
  addObsColumn(colSchema, Ctor, value) {
    // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
    const annoMatrix = (this as any).annoMatrix.addObsColumn(
      colSchema,
      Ctor,
      value
    );
    // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
    const obsCrossfilter = (this as any).obsCrossfilter.setData(
      annoMatrix._cache.obs
    );
    return new AnnoMatrixObsCrossfilter(annoMatrix, obsCrossfilter);
  }

  // @ts-expect-error ts-migrate(7006) FIXME: Parameter 'col' implicitly has an 'any' type.
  // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types --- FIXME: disabled temporarily on migrate to TS.
  dropObsColumn(col) {
    // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
    const annoMatrix = (this as any).annoMatrix.dropObsColumn(col);
    // @ts-expect-error ts-migrate(2339) FIXME: Property 'obsCrossfilter' does not exist on type '... Remove this comment to see the full error message
    let { obsCrossfilter } = this;
    const dimName = _dimensionName("obs", col);
    if (obsCrossfilter.hasDimension(dimName)) {
      obsCrossfilter = obsCrossfilter.delDimension(dimName);
    }
    return new AnnoMatrixObsCrossfilter(annoMatrix, obsCrossfilter);
  }

  // @ts-expect-error ts-migrate(7006) FIXME: Parameter 'oldCol' implicitly has an 'any' type.
  // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types --- FIXME: disabled temporarily on migrate to TS.
  renameObsColumn(oldCol, newCol) {
    // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
    const annoMatrix = (this as any).annoMatrix.renameObsColumn(oldCol, newCol);
    const oldDimName = _dimensionName("obs", oldCol);
    const newDimName = _dimensionName("obs", newCol);
    // @ts-expect-error ts-migrate(2339) FIXME: Property 'obsCrossfilter' does not exist on type '... Remove this comment to see the full error message
    let { obsCrossfilter } = this;
    if (obsCrossfilter.hasDimension(oldDimName)) {
      obsCrossfilter = obsCrossfilter.renameDimension(oldDimName, newDimName);
    }
    return new AnnoMatrixObsCrossfilter(annoMatrix, obsCrossfilter);
  }

  // @ts-expect-error ts-migrate(7006) FIXME: Parameter 'col' implicitly has an 'any' type.
  // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types --- FIXME: disabled temporarily on migrate to TS.
  addObsAnnoCategory(col, category) {
    // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
    const annoMatrix = (this as any).annoMatrix.addObsAnnoCategory(
      col,
      category
    );
    const dimName = _dimensionName("obs", col);
    // @ts-expect-error ts-migrate(2339) FIXME: Property 'obsCrossfilter' does not exist on type '... Remove this comment to see the full error message
    let { obsCrossfilter } = this;
    if (obsCrossfilter.hasDimension(dimName)) {
      obsCrossfilter = obsCrossfilter.delDimension(dimName);
    }
    return new AnnoMatrixObsCrossfilter(annoMatrix, obsCrossfilter);
  }

  // @ts-expect-error ts-migrate(7006) FIXME: Parameter 'col' implicitly has an 'any' type.
  // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types --- FIXME: disabled temporarily on migrate to TS.
  async removeObsAnnoCategory(col, category, unassignedCategory) {
    // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
    const annoMatrix = await (this as any).annoMatrix.removeObsAnnoCategory(
      col,
      category,
      unassignedCategory
    );
    const dimName = _dimensionName("obs", col);
    // @ts-expect-error ts-migrate(2339) FIXME: Property 'obsCrossfilter' does not exist on type '... Remove this comment to see the full error message
    let { obsCrossfilter } = this;
    if (obsCrossfilter.hasDimension(dimName)) {
      obsCrossfilter = obsCrossfilter.delDimension(dimName);
    }
    return new AnnoMatrixObsCrossfilter(annoMatrix, obsCrossfilter);
  }

  // @ts-expect-error ts-migrate(7006) FIXME: Parameter 'col' implicitly has an 'any' type.
  // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types --- FIXME: disabled temporarily on migrate to TS.
  async setObsColumnValues(col, rowLabels, value) {
    // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
    const annoMatrix = await (this as any).annoMatrix.setObsColumnValues(
      col,
      rowLabels,
      value
    );
    const dimName = _dimensionName("obs", col);
    // @ts-expect-error ts-migrate(2339) FIXME: Property 'obsCrossfilter' does not exist on type '... Remove this comment to see the full error message
    let { obsCrossfilter } = this;
    if (obsCrossfilter.hasDimension(dimName)) {
      obsCrossfilter = obsCrossfilter.delDimension(dimName);
    }
    return new AnnoMatrixObsCrossfilter(annoMatrix, obsCrossfilter);
  }

  // @ts-expect-error ts-migrate(7006) FIXME: Parameter 'col' implicitly has an 'any' type.
  // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types --- FIXME: disabled temporarily on migrate to TS.
  async resetObsColumnValues(col, oldValue, newValue) {
    // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
    const annoMatrix = await (this as any).annoMatrix.resetObsColumnValues(
      col,
      oldValue,
      newValue
    );
    const dimName = _dimensionName("obs", col);
    // @ts-expect-error ts-migrate(2339) FIXME: Property 'obsCrossfilter' does not exist on type '... Remove this comment to see the full error message
    let { obsCrossfilter } = this;
    if (obsCrossfilter.hasDimension(dimName)) {
      obsCrossfilter = obsCrossfilter.delDimension(dimName);
    }
    return new AnnoMatrixObsCrossfilter(annoMatrix, obsCrossfilter);
  }

  // @ts-expect-error ts-migrate(7006) FIXME: Parameter 'colSchema' implicitly has an 'any' type... Remove this comment to see the full error message
  // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types --- FIXME: disabled temporarily on migrate to TS.
  addEmbedding(colSchema) {
    // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
    const annoMatrix = (this as any).annoMatrix.addEmbedding(colSchema);
    return new AnnoMatrixObsCrossfilter(
      annoMatrix,
      // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
      (this as any).obsCrossfilter
    );
  }

  /**
   * Drop the crossfilter dimension. Do not change the annoMatrix. Useful when we
   * want to stop trackin the selection state, but aren't sure we want to blow the
   * annomatrix cache.
   */
  // @ts-expect-error ts-migrate(7006) FIXME: Parameter 'field' implicitly has an 'any' type.
  // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types --- FIXME: disabled temporarily on migrate to TS.
  dropDimension(field, query) {
    // @ts-expect-error ts-migrate(2339) FIXME: Property 'annoMatrix' does not exist on type 'Anno... Remove this comment to see the full error message
    const { annoMatrix } = this;
    // @ts-expect-error ts-migrate(2339) FIXME: Property 'obsCrossfilter' does not exist on type '... Remove this comment to see the full error message
    let { obsCrossfilter } = this;
    const keys = annoMatrix
      .getCacheKeys(field, query)
      // @ts-expect-error ts-migrate(7006) FIXME: Parameter 'k' implicitly has an 'any' type.
      .filter((k) => k !== undefined);
    const dimName = _dimensionName(field, keys);
    if (obsCrossfilter.hasDimension(dimName)) {
      obsCrossfilter = obsCrossfilter.delDimension(dimName);
    }
    return new AnnoMatrixObsCrossfilter(annoMatrix, obsCrossfilter);
  }

  /**
  Selection state - API is identical to ImmutableTypedCrossfilter, as these
  are just wrappers to lazy create indices.
  **/

  // @ts-expect-error ts-migrate(7006) FIXME: Parameter 'field' implicitly has an 'any' type.
  // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types --- FIXME: disabled temporarily on migrate to TS.
  async select(field, query, spec) {
    // @ts-expect-error ts-migrate(2339) FIXME: Property 'annoMatrix' does not exist on type 'Anno... Remove this comment to see the full error message
    const { annoMatrix } = this;
    // @ts-expect-error ts-migrate(2339) FIXME: Property 'obsCrossfilter' does not exist on type '... Remove this comment to see the full error message
    let { obsCrossfilter } = this;

    if (!annoMatrix?._cache?.[field]) {
      throw new Error("Unknown field name");
    }
    if (field === "var") {
      throw new Error("unable to obsSelect upon the var dimension");
    }

    // grab the data, so we can grab the index.
    const df = await annoMatrix.fetch(field, query);

    const dimName = _dimensionNameFromDf(field, df);
    if (!obsCrossfilter.hasDimension(dimName)) {
      // lazy index generation - add dimension when first used
      obsCrossfilter = this._addObsCrossfilterDimension(
        annoMatrix,
        obsCrossfilter,
        field,
        df
      );
    }

    // select
    obsCrossfilter = obsCrossfilter.select(dimName, spec);
    return new AnnoMatrixObsCrossfilter(annoMatrix, obsCrossfilter);
  }

  // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types --- FIXME: disabled temporarily on migrate to TS.
  selectAll() {
    /*
		Select all on any dimension in this field.
		*/
    // @ts-expect-error ts-migrate(2339) FIXME: Property 'annoMatrix' does not exist on type 'Anno... Remove this comment to see the full error message
    const { annoMatrix } = this;
    // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
    const currentDims = (this as any).obsCrossfilter.dimensionNames();
    const obsCrossfilter = currentDims.reduce(
      // @ts-expect-error ts-migrate(7006) FIXME: Parameter 'xfltr' implicitly has an 'any' type.
      (xfltr, dim) => xfltr.select(dim, { mode: "all" }),
      // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
      (this as any).obsCrossfilter
    ); // eslint-disable-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
    return new AnnoMatrixObsCrossfilter(annoMatrix, obsCrossfilter);
  }

  // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types --- FIXME: disabled temporarily on migrate to TS.
  countSelected() {
    /* if no data yet indexed in the crossfilter, just say everything is selected */
    // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
    if ((this as any).obsCrossfilter.size() === 0)
      // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
      return (this as any).annoMatrix.nObs;
    // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
    return (this as any).obsCrossfilter.countSelected();
  }

  // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types --- FIXME: disabled temporarily on migrate to TS.
  allSelectedMask() {
    /* if no data yet indexed in the crossfilter, just say everything is selected */
    if (
      // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
      (this as any).obsCrossfilter.size() === 0 ||
      // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
      (this as any).obsCrossfilter.dimensionNames().length === 0
    ) {
      /* fake the mask */
      // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
      return new Uint8Array((this as any).annoMatrix.nObs).fill(1);
    }
    // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
    return (this as any).obsCrossfilter.allSelectedMask();
  }

  // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types --- FIXME: disabled temporarily on migrate to TS.
  allSelectedLabels() {
    /* if no data yet indexed in the crossfilter, just say everything is selected */
    if (
      // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
      (this as any).obsCrossfilter.size() === 0 ||
      // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
      (this as any).obsCrossfilter.dimensionNames().length === 0
    ) {
      // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
      return (this as any).annoMatrix.rowIndex.labels();
    }

    // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
    const mask = (this as any).obsCrossfilter.allSelectedMask();
    // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
    const index = (this as any).annoMatrix.rowIndex.isubsetMask(mask);
    return index.labels();
  }

  // @ts-expect-error ts-migrate(7006) FIXME: Parameter 'array' implicitly has an 'any' type.
  // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types --- FIXME: disabled temporarily on migrate to TS.
  fillByIsSelected(array, selectedValue, deselectedValue) {
    /* if no data yet indexed in the crossfilter, just say everything is selected */
    if (
      // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
      (this as any).obsCrossfilter.size() === 0 ||
      // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
      (this as any).obsCrossfilter.dimensionNames().length === 0
    ) {
      return array.fill(selectedValue);
    }
    // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
    return (this as any).obsCrossfilter.fillByIsSelected(
      array,
      selectedValue,
      deselectedValue
    );
  }

  /**
   ** Private below
   **/

  // @ts-expect-error ts-migrate(7006) FIXME: Parameter 'annoMatrix' implicitly has an 'any' typ... Remove this comment to see the full error message
  // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types --- FIXME: disabled temporarily on migrate to TS.
  _addObsCrossfilterDimension(annoMatrix, obsCrossfilter, field, df) {
    if (field === "var") return obsCrossfilter;
    const dimName = _dimensionNameFromDf(field, df);
    const dimParams = this._getObsDimensionParams(field, df);
    obsCrossfilter = obsCrossfilter.setData(annoMatrix._cache.obs);
    // @ts-expect-error ts-migrate(2488) FIXME: Type 'any[] | undefined' must have a '[Symbol.iter... Remove this comment to see the full error message
    obsCrossfilter = obsCrossfilter.addDimension(dimName, ...dimParams);
    return obsCrossfilter;
  }

  // @ts-expect-error ts-migrate(7006) FIXME: Parameter 'field' implicitly has an 'any' type.
  // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types --- FIXME: disabled temporarily on migrate to TS.
  _getColumnBaseType(field, col) {
    /* Look up the primitive type for this field/col */
    const colSchema = _getColumnSchema(
      // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
      (this as any).annoMatrix.schema,
      field,
      col
    );
    return colSchema.type;
  }

  // @ts-expect-error ts-migrate(7006) FIXME: Parameter 'field' implicitly has an 'any' type.
  // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types --- FIXME: disabled temporarily on migrate to TS.
  _getObsDimensionParams(field, df) {
    /* return the crossfilter dimensiontype type and params for this field/dataframe */

    if (field === "emb") {
      /* assumed to be 2D */
      return ["spatial", df.icol(0).asArray(), df.icol(1).asArray()];
    }

    /* assumed to be 1D */
    const col = df.icol(0);
    const colName = df.colIndex.getLabel(0);
    const type = this._getColumnBaseType(field, colName);
    if (type === "string" || type === "categorical" || type === "boolean") {
      return ["enum", col.asArray()];
    }
    if (type === "int32") {
      return ["scalar", col.asArray(), Int32Array];
    }
    if (type === "float32") {
      return ["scalar", col.asArray(), Float32Array];
    }
    // Currently not supporting boolean and categorical types.
    console.error(
      `Warning - unknown metadata schema (${type}) for field ${field} ${colName}.`
    );
    // skip it - we don't know what to do with this type

    return undefined;
  }
}
