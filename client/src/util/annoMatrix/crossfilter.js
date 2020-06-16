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

  annoMatrix() {
    return this.annoMatrix();
  }

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

    const dimName = dimensionName(field, df);
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
    const dimName = dimensionName(field, df);
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

function dimensionName(field, df) {
  const colNames = df.colIndex.labels().join(":");
  return `${field}/${colNames}`;
}
