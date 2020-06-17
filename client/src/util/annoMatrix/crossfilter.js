/*
Row crossfilter proxy for an AnnoMatrix.  This wraps Crossfilter,
providing a number of services:
	- on-demand index creation as data is loaded
	- transparently mapping between queries and crossfilter index names
*/
import Crossfilter from "../typedCrossfilter";
import { getColumnSchema } from "./schema";

export default class AnnoMatrixObsCrossfilter {
  constructor(annoMatrix, obsCrossfilter = null) {
    this.annoMatrix = annoMatrix;
    this.obsCrossfilter = obsCrossfilter || new Crossfilter(annoMatrix.obs);
  }

  size() {
    return this.obsCrossfilter.size();
  }

  /**
  Managing the associated annoMatrix.  These wrappers are necessary to 
  make coordinated changes to BOTH the crossfilter and annoMatrix, and
  ensure that all state stays synchronized.
  **/
  addObsColumn(colSchema, Ctor, value) {
    const annoMatrix = this.annoMatrix.addObsColumn(colSchema, Ctor, value);
    const obsCrossfilter = this.obsCrossfilter.setData(annoMatrix.obs);
    return new AnnoMatrixObsCrossfilter(annoMatrix, obsCrossfilter);
  }

  dropObsColumn(col) {
    const annoMatrix = this.annoMatrix.dropObsColumn(col);
    let { obsCrossfilter } = this;
    const dimName = dimensionName("obs", col);
    if (obsCrossfilter.hasDimension(dimName)) {
      obsCrossfilter = obsCrossfilter.delDimension(dimName);
    } else {
      obsCrossfilter = obsCrossfilter.setData(annoMatrix.obs);
    }
    return new AnnoMatrixObsCrossfilter(annoMatrix, obsCrossfilter);
  }

  renameObsColumn(oldCol, newCol) {
    const annoMatrix = this.annoMatrix.renameObsColumn(oldCol, newCol);
    const oldDimName = dimensionName("obs", oldCol);
    const newDimName = dimensionName("obs", newCol);
    let { obsCrossfilter } = this;
    if (obsCrossfilter.hasDimension(oldDimName)) {
      obsCrossfilter = obsCrossfilter.renameDimension(oldDimName, newDimName);
    } else {
      obsCrossfilter = obsCrossfilter.setData(annoMatrix.obs);
    }
    return new AnnoMatrixObsCrossfilter(annoMatrix, obsCrossfilter);
  }

  addObsAnnoCategory(col, category) {
    const annoMatrix = this.annoMatrix.addObsAnnoCategory(col, category);
    const dimName = dimensionName("obs", col);
    let { obsCrossfilter } = this;
    if (obsCrossfilter.hasDimension(dimName)) {
      obsCrossfilter = obsCrossfilter.delDimension(dimName);
    }
    return new AnnoMatrixObsCrossfilter(annoMatrix, obsCrossfilter);
  }

  removeObsAnnoCategory(col, category) {
    const annoMatrix = this.annoMatrix.removeObsAnnoCategory(col, category);
    const dimName = dimensionName("obs", col);
    let { obsCrossfilter } = this;
    if (obsCrossfilter.hasDimension(dimName)) {
      obsCrossfilter = obsCrossfilter.delDimension(dimName);
    }
    return new AnnoMatrixObsCrossfilter(annoMatrix, obsCrossfilter);
  }

  setObsColumnValues(col, rowLabels, value) {
    const annoMatrix = this.annoMatrix.setObsColumnValues(
      col,
      rowLabels,
      value
    );
    const dimName = dimensionName("obs", col);
    let { obsCrossfilter } = this;
    if (obsCrossfilter.hasDimension(dimName)) {
      obsCrossfilter = obsCrossfilter.delDimension(dimName);
    }
    return new AnnoMatrixObsCrossfilter(annoMatrix, obsCrossfilter);
  }
  /**
  Selection state 
  **/

  async select(field, query, spec) {
    const { annoMatrix } = this;
    let { obsCrossfilter } = this;

    if (!annoMatrix?.[field]) {
      throw new Error("Unknown field name");
    }
    if (field === "var") {
      throw new Error("unable to obsSelect upon the var dimension");
    }

    // grab the data, so we can grab the index.
    const df = await annoMatrix.fetch(field, query);

    const dimName = dimensionNameFromDf(field, df);
    if (!obsCrossfilter.hasDimension(dimName)) {
      // lazy index generation - add dimension when first used
      obsCrossfilter = this._addObsCrossfilterDimension(
        obsCrossfilter,
        field,
        df
      );
    }

    // select
    obsCrossfilter = obsCrossfilter.select(dimName, spec);
    return new AnnoMatrixObsCrossfilter(annoMatrix, obsCrossfilter);
  }

  selectAll() {
    /*
		Select all on any dimension in this field.
		*/
    const { annoMatrix } = this;
    const currentDims = this.obsCrossfilter.dimensionNames();
    const obsCrossfilter = currentDims.reduce((xfltr, dim) => {
      return xfltr.select(dim, { mode: "all" });
    }, this.obsCrossfilter);
    return new AnnoMatrixObsCrossfilter(annoMatrix, obsCrossfilter);
  }

  countSelected() {
    /* if no data yet indexed in the crossfilter, just say everything is selected */
    if (this.obsCrossfilter.size() === 0) return this.annoMatrix.nObs;
    return this.obsCrossfilter.countSelected();
  }

  allSelectedMask() {
    /* if no data yet indexed in the crossfilter, just say everything is selected */
    if (
      this.obsCrossfilter.size() === 0 ||
      this.obsCrossfilter.dimensionNames().length === 0
    ) {
      /* fake the mask */
      return new Uint8Array(this.annoMatrix.nObs).fill(1);
    }
    return this.obsCrossfilter.allSelectedMask();
  }

  allSelectedLabels() {
    /* if no data yet indexed in the crossfilter, just say everything is selected */
    if (
      this.obsCrossfilter.size() === 0 ||
      this.obsCrossfilter.dimensionNames().length === 0
    ) {
      return this.annoMatrix.rowIndex;
    }

    const mask = this.obsCrossfilter.allSelectedMask();
    const index = this.annoMatrix.rowIndex.isubsetMask(mask);
    return index.labels();
  }

  fillByIsSelected(array, selectedValue, deselectedValue) {
    /* if no data yet indexed in the crossfilter, just say everything is selected */
    if (
      this.obsCrossfilter.size() === 0 ||
      this.obsCrossfilter.dimensionNames().length === 0
    ) {
      return new Uint8Array(this.annoMatrix.nObs).fill(selectedValue);
    }
    return this.obsCrossfilter.fillByIsSelected(
      array,
      selectedValue,
      deselectedValue
    );
  }

  /**
   ** Private
   **/

  _addObsCrossfilterDimension(obsCrossfilter, field, df) {
    if (field === "var") return obsCrossfilter;
    const dimName = dimensionNameFromDf(field, df);
    const dimParams = this._getObsDimensionParams(field, df);
    if (obsCrossfilter.size() === 0)
      obsCrossfilter = obsCrossfilter.setData(df);
    obsCrossfilter = obsCrossfilter.addDimension(dimName, ...dimParams);
    return obsCrossfilter;
  }

  _getColumnBaseType(field, col) {
    /* Look up the primitive type for this field/col */
    const colSchema = getColumnSchema(this.annoMatrix.schema, field, col);
    return colSchema.type;
  }

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

function dimensionNameFromDf(field, df) {
  const colNames = df.colIndex.labels();
  return dimensionName(field, colNames);
}

function dimensionName(field, colNames) {
  if (!Array.isArray(colNames)) return `${field}/${colNames}`;
  return `${field}/${colNames.join(":")}`;
}
