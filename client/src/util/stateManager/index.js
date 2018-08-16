"use strict";
// jshint esversion: 6

import _ from "lodash";

/*
Model manager providing an abstraction for the use of the reducer code.
This module provides several buckets of functionality:
  - schema and config driven tranformation of the dataframe wire protocol
    into a format that is easy for the UI code to use.
  - manage the universe/world abstraction:
    + universe: all of the server-provided, read-only data
    + world: subset of universe
  - lazy access and caching of dataframe contents as needed
*/

class Universe {
  /*
  creates empty universe.
  alternative: should the constructor consume the raw data needed
  to create a universe, removing the need for the setters below?
  */
  constructor(apiVersion) {}

  /**
   ** setters or equivalent - build the universe
   **
   ** to be used by master reducer
   **/

  /*
  Init-time set server-provided configuration, eg,
  setSchema(action.data).
  */
  setSchema(OTAresponse) {}
  setObsAnnotationsObs(OTAresponse) {}
  setAnnotationsVar(OTAresponse) {}
  setLayoutObs(OTAresponse) {}
  // more?

  /*
  Universe is built - freeze and generate any derivative information we
  compute on the client side (eg, ranges).
  */
  finalize() {}

  /**
   ** gettors/accessors
   **
   ** Universe doesn't get used directly by UI code - it should use world.
   **
   ** some of these getters may be replaced by actual objects, rather than
   ** getter functions - listing them as getters to make interface clear.
   **/

  /* return schema in same format as REST 0.2 API */
  get schema() {}
}

/*
World is a subset of universe.   All application UI should use world, and should
(generally) not use Universe.
*/
class World {
  /**
   ** World "creators" - these each create a world, but in different ways
   **/

  /*
  create world which is eq universe.  Universe MUST be fully
  initialized for this to suceed (eg all setters called, finalize()
  called)
  */
  constructor(universe) {}

  /*
  create a world that is a subset of first argument, where the new world contains
  a subset of universe defined by 'createFrom' argument.

  'createFrom' argument can be one of:
    * 'universe': create world == universe.  Alias for new World(universe)
    * 'selected': create world containing obs/cells currently selected in
      this world.crossfilter
  Others may be added in the future.
  */
  contructor(world, createFrom) {}

  /**
   ** Getters for the state inside world.    These are guarnateed to return
   ** data for just the subset included in world.
   **/

  /*
  get schema.  Literally an alias for universe.schema as schema doesn't change
  when world changes.
  */
  get schema() {}

  /* the crossfilter object for this world */
  get crossfilter() {}

  /*  dimensions for annotation/metadata, orgaized as an object { 'fieldName': dimension, ... } */
  get dimensionMap() {}

  /*
  annotations, layout  and matrix getters

  annotations will return an array of objects.  Each object contains all annotation
  values for a given observation/cell, keyed by annotation name, PLUS a key
  '__cellId__', containing a REST API ID for this obs/cell (referred to as the
  obsIndex in the REST 0.2 spec or cellIndex in the 0.1 spec.

  Example:  [ { __cellId__: 99, cluster: 'blue', numReads: 93933 } ]

  NOTE: world.annotation should be identical to the old state.cells value,
  EXCEPT that
    * __cellIndex__ renamed to __cellId__  XXX: should be in mapToUniverse
    * __x__ and __y__ are now in world.layout
    * __color__ and __colorRBG__ should be moved to controls reducer


  layout will return an object containing two arrays, containing X and Y
  coordinates respectively.

  Example: { X: [ 0.33, 0.23, ... ], Y: [ 0.8, 0.777, ... ]}


  matrix:  TBD


  Array position will be consistent across annotations and layout return vals,
  so world.annotations[0] and world.layout.X[0] refer to the same obs/cell.

  */
  get annotations() {}
  get layout() {}

  // private state - maps back to Universe index for any given obs
  get mapToUniverseIndex() {}

  /*
  Summary information for each obs/cell annotation, keyed by annotation name.
  Value will be an object, containing either 'range' or 'options' object,
  depending on the annotation schema type (categorical or continuous)

  XXX: code needs a list of all gene names, eg, "names": ...

  Example:
    {
      "Splice_sites_Annotated": {
        "range": {
          "min": 26,
          "max": 1075869
        }
      },
      "Selection": {
        "options": {
          "Astrocytes(HEPACAM)": 714,
          "Endothelial(BSC)": 123,
          "Oligodendrocytes(GC)": 294,
          "Neurons(Thy1)": 685,
          "Microglia(CD45)": 1108,
          "Unpanned": 665
        }
      }
    }
  */
  get summary() {}
}

/*
This is how I imagine the reducer will use this code.
*/

const DataFrame = (
  state = {
    universe: null,
    world: null,
    initState: { schema: false, annotations: false, layout: false }
  },
  action
) => {
  let isInitAction = false;
  const initState = state.initState;
  const universe = state.universe ? state.universe : new Universe();

  /*
  First, handle initialization related actions
  */
  switch (action.type) {
    case "initialize - schema fetch success": {
      universe.setSchemaFromOTA(action.data);
      initState.schema = true;
      isInitAction = true;
      break;
    }
    case "initialize - all obs annotation fetch success": {
      universe.setAnnotationsFromOTA(action.data);
      initState.annotations = true;
      isInitAction = true;
      break;
    }
    case "initialize - all obs layout fetch success": {
      universe.setLayoutFromOTA(action.data);
      initState.layout = true;
      isInitAction = true;
      break;
    }
  }
  if (isInitAction && _.every(state.initState)) {
    universe.finalize();
    return Object.assign({}, state, {
      universe,
      world: new World(universe),
      initState
    });
  }

  /*
  Second, handle changes to world definition, which only work once we are init
  */
  switch (action.type) {
    case "reset world to be all of the universe": {
      // this action currently called "reset graph"
      return Object.assign({}, state, {
        world: new World(state.world, "universe")
      });
    }
    case "set world to be current selection": {
      return Object.assign({}, state, {
        world: new World(state.world, "selection")
      });
    }
    /* NOTE: we can add other modes to the World constructor that
    support creating a new World from other "collections" of cells. */
  }
};
