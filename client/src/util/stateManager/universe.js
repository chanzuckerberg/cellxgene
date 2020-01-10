import _ from "lodash";

import { unassignedCategoryLabel } from "../../globals";
import { decodeMatrixFBS } from "./matrix";
import * as Dataframe from "../dataframe";
import { isFpTypedArray } from "../typeHelpers";
import { indexEntireSchema, sortAllCategorical } from "./schemaHelpers";
import { isCategoricalAnnotation } from "./annotationsHelpers";

/*
Private helper function - create and return a template Universe
*/
function templateUniverse() {
  /* default universe template */
  return {
    nObs: 0,
    nVar: 0,
    schema: {},

    /*
    annotations
    */
    obsAnnotations: Dataframe.Dataframe.empty(),
    varAnnotations: Dataframe.Dataframe.empty(),
    /*
    layout
    */
    obsLayout: Dataframe.Dataframe.empty(),

    /*
    Var data columns - subset of all
    */
    varData: Dataframe.Dataframe.empty(null, new Dataframe.KeyIndex())
  };
}

/*
This module implements functions that support storage of "Universe",
aka all of the var/obs data and annotations.

These functions are used exclusively by the actions and reducers to
build an internal POJO for use by the rendering components.
*/

function promoteTypedArray(o) {
  /*
  Decide what internal data type to use for the data returned from 
  the server.

  TODO - future optimization: not all int32/uint32 data series require
  promotion to float64.  We COULD simply look at the data to decide. 
  */
  if (isFpTypedArray(o) || Array.isArray(o)) return o;

  let TyepdArrayCtor;
  switch (o.constructor) {
    case Int8Array:
    case Uint8Array:
    case Uint8ClampedArray:
    case Int16Array:
    case Uint16Array:
      TyepdArrayCtor = Float32Array;
      break;

    case Int32Array:
    case Uint32Array:
      TyepdArrayCtor = Float64Array;
      break;

    default:
      throw new Error("Unexpected data type returned from server.");
  }
  if (o.constructor === TyepdArrayCtor) return o;
  return new TyepdArrayCtor(o);
}

function AnnotationsFBSToDataframe(arrayBuffers) {
  /*
  Convert a Matrix FBS to a Dataframe.

  The application has strong assumptions that all scalar data will be
  stored as a float32 or float64 (regardless of underlying data types).
  For example, clipping of value ranges (eg, user-selected percentiles)
  depends on the ability to use NaN in any numeric type.

  All float data from the server is left as is.  All non-float is promoted
  to an appropriate float.
  */
  if (arrayBuffers instanceof ArrayBuffer) {
    arrayBuffers = [arrayBuffers];
  }
  if (arrayBuffers.length === 0) {
    return Dataframe.Dataframe.empty();
  }

  const fbs = arrayBuffers.map(ab => decodeMatrixFBS(ab, true)); // leave in place
  /* check that all FBS have same row dimensionality */
  const nRows = fbs[0].nRows;
  fbs.forEach(b => {
    if (b.nRows !== nRows)
      throw new Error("FBS with inconsistent dimensionality");
  });
  const columns = fbs
    .map(fb =>
      fb.columns.map(c => {
        if (isFpTypedArray(c) || Array.isArray(c)) return c;
        return promoteTypedArray(c);
      })
    )
    .flat();
  const colIdx = fbs.map(b => b.colIdx).flat();
  const nCols = columns.length;

  const df = new Dataframe.Dataframe(
    [nRows, nCols],
    columns,
    null,
    new Dataframe.KeyIndex(colIdx)
  );
  return df;
}

function reconcileSchemaCategoriesWithSummary(universe) {
  /*
  where we treat types as (essentially) categorical metadata, update
  the schema with data-derived categories (in addition to those in
  the server declared schema).

  For example, boolean defined fields in the schema do not contain
  explicit declaration of categories (nor do string fields).  In these
  cases, add a 'categories' field to the schema so it is accessible.

  In addition, we have a client-side convention (UI) that all writable
  annotations must have an 'unassigned' category, even if it is not currently
  in use.
  */

  universe.schema.annotations.obs.columns.forEach(s => {
    const { type, name, writable } = s;
    if (type === "string" || type === "boolean" || type === "categorical") {
      if (universe.obsAnnotations.hasCol(name)) {
        const categories = _.union(
          s.categories ?? [],
          universe.obsAnnotations.col(name).summarize().categories ?? []
        );
        s.categories = categories;
      }
    }

    if (writable && s.categories.indexOf(unassignedCategoryLabel) === -1) {
      s.categories = s.categories.concat(unassignedCategoryLabel);
    }
  });
}

export function createUniverseFromResponse(
  configResponse,
  schemaResponse,
  annotationsObsResponse,
  annotationsVarResponse,
  layoutFBSResponse
) {
  /*
  build & return universe from a REST 0.2 /config, /schema and /annotations/obs response
  */
  const { schema } = schemaResponse;
  const universe = templateUniverse();

  /* schema related */
  universe.schema = schema;
  universe.nObs = schema.dataframe.nObs;
  universe.nVar = schema.dataframe.nVar;
  /* add defaults, as we can't assume back-end will fully populate schema */
  if (!schema.layout.var) schema.layout.var = [];
  if (!schema.layout.obs) schema.layout.obs = [];

  /* annotations */
  universe.obsAnnotations = AnnotationsFBSToDataframe(annotationsObsResponse);
  universe.varAnnotations = AnnotationsFBSToDataframe(annotationsVarResponse);
  /* layout */
  // universe.obsLayout = LayoutFBSToDataframe(layoutFBSResponse);
  universe.obsLayout = AnnotationsFBSToDataframe(layoutFBSResponse);

  /* sanity checks */
  if (
    universe.nObs !== universe.obsLayout.length ||
    universe.nObs !== universe.obsAnnotations.length ||
    universe.nVar !== universe.varAnnotations.length
  ) {
    throw new Error("Universe dimensionality mismatch - failed to load");
  }

  reconcileSchemaCategoriesWithSummary(universe);
  sortAllCategorical(universe.schema);
  indexEntireSchema(universe.schema);

  /* sanity checks */
  if (
    schema.annotations.obs.columns.some(
      s => s.writable && !isCategoricalAnnotation(schema, s.name)
    )
  ) {
    throw new Error(
      "Writable continuous obs annotations are not supported - failed to load"
    );
  }

  return universe;
}

export function convertDataFBStoObject(universe, arrayBuffer) {
  /*
  /data/var returns a flatbuffer (FBS) as described by cellxgene/fbs/matrix.fbs

  This routine converts the binary wire encoding into a JS object:

  {
    gene: Float32Array,
    ...
  }
  */
  const fbs = decodeMatrixFBS(arrayBuffer);
  const { colIdx, columns } = fbs;
  const result = {};

  if (!columns.every(isFpTypedArray)) {
    // We have strong assumptions that all var data is float
    throw new Error("Unexpected non-floating point response from server.");
  }

  const varIndexName = universe.schema.annotations.var.index;
  for (let c = 0; c < colIdx.length; c += 1) {
    const varName = universe.varAnnotations.at(colIdx[c], varIndexName);
    result[varName] = columns[c];
  }
  return result;
}
