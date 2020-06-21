// jshint esversion: 6

/*
Model manager providing an abstraction for the use of the reducer code.
This module provides several buckets of functionality:
  - schema and config driven tranformation of the wire protocol
    into a format that is easy for the UI code to use.
  - manage the universe/world abstraction:
    + universe: all of the server-provided, read-only data
    + world: subset of universe
  - lazy access and caching of dataframe contents as needed

This is all VERY tightly integrated with reducers and actions, and
exists to support those concepts.
*/

export * as ColorHelpers from "./colorHelpers";
export * as ControlsHelpers from "./controlsHelpers";
export * as AnnotationsHelpers from "./annotationsHelpers";
export * as SchemaHelpers from "./schemaHelpers";
export * as MatrixFBS from "./matrix";
