// jshint esversion: 6

import _ from "lodash";

import summarizeAnnotations from "./summarizeAnnotations";
import decodeMatrixFBS from "./matrix";
import * as Dataframe from "../dataframe";

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
    Annotations
    */
    obsAnnotations: Dataframe.Dataframe.empty(),
    varAnnotations: Dataframe.Dataframe.empty(),
    obsLayout: Dataframe.Dataframe.empty(),
    summary: null /* derived data summaries. XXX: consider exploding in place */,

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

function AnnotationsFBSToDataframe(arrayBuffer) {
  /*
  Convert a Matrix FBS to a Dataframe.
  */
  const fbs = decodeMatrixFBS(arrayBuffer);
  const df = new Dataframe.Dataframe(
    [fbs.nRows, fbs.nCols],
    fbs.columns,
    null,
    new Dataframe.KeyIndex(fbs.colIdx)
  );
  return df;
}

function LayoutFBSToDataframe(arrayBuffer) {
  const fbs = decodeMatrixFBS(arrayBuffer, true);
  const df = new Dataframe.Dataframe(
    [fbs.nRows, fbs.nCols],
    fbs.columns,
    null,
    new Dataframe.KeyIndex(["X", "Y"])
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
  */

  _.forEach(universe.schema.annotations.obs, s => {
    if (
      s.type === "string" ||
      s.type === "boolean" ||
      s.type === "categorical"
    ) {
      const categories = _.union(
        _.get(s, "categories", []),
        _.get(universe.summary.obs[s.name], "categories", [])
      );
      s.categories = categories;
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

  /* annotations */
  universe.obsAnnotations = AnnotationsFBSToDataframe(annotationsObsResponse);
  universe.varAnnotations = AnnotationsFBSToDataframe(annotationsVarResponse);
  /* layout */
  universe.obsLayout = LayoutFBSToDataframe(layoutFBSResponse);

  /* sanity check */
  if (
    universe.nObs !== universe.obsLayout.length ||
    universe.nObs !== universe.obsAnnotations.length ||
    universe.nVar !== universe.varAnnotations.length
  ) {
    throw new Error("Universe dimensionality mismatch - failed to load");
  }

  universe.summary = summarizeAnnotations(
    universe.schema,
    universe.obsAnnotations,
    universe.varAnnotations
  );

  reconcileSchemaCategoriesWithSummary(universe);
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

  for (let c = 0; c < colIdx.length; c += 1) {
    const varName = universe.varAnnotations.at(colIdx[c], "name");
    result[varName] = columns[c];
  }
  return result;
}
