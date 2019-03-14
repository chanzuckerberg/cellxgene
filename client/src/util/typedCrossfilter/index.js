/*
Crossfilter - a re-implementation of a subset of crossfilter, with
time/space optimizations predicated upon the following assumptions:
  - dimensions are uniformly typed, and all values must be of that type
  - dimension values must be a primitive type (int, float, string).  Arrays
    or other complex types not supported.
  - dimension creation requires call-provided type declaration
  - no support for adding/removing data to an existing crossfilter.  If you
    want to do that, you have to create the new crossfilter, using the new
    data, from scratch.

In addition, this implementation is easier to use with a "redux" style
app, as all operations on the crossfilter are immutable (ie, return a
new crossfilter).

The actual backing store for a dimension is a TypedArray, enabling
performance improvements over the original crossfilter.

There are also a handful of new methods, primarily to take advantage of the
performance (eg, crossfilter.fillBySelection)

Helpful documents (this module follows similar concepts as the original,
but deviates from the API):
    https://github.com/square/crossfilter/
    http://square.github.io/crossfilter/

There is also a newer, community supported fork of crossfilter, with a
more complex API.  In a few cases, elements of that API were incorporated.
    https://github.com/square/crossfilter/

See test cases for some concrete examples.
*/

export { default } from "./crossfilter";
