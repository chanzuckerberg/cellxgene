import { IdentityInt32Index, isLabelIndex } from "./labelIndex";
// weird cross-dependency that we should clean up someday...
import { callOnceLazy, memoize, __getMemoId } from "./util";
import { isTypedArray, isAnyArray } from "../../common/types/arraytypes";
import {
  summarizeContinuous,
  summarizeCategorical as _summarizeCategorical,
} from "./summarize";
import {
  histogramCategorical as _histogramCategorical,
  hashCategorical,
  histogramContinuous,
  hashContinuous,
} from "./histogram";

/*
Dataframe is an immutable 2D matrix similiar to Python Pandas Dataframe,
but (currently) without all of the surrounding support functions.
Data is stored in column-major layout, and each column is monomorphic.

It supports:
* Relatively efficient creation, cloning and subsetting
* Very efficient columnar access (eg, sum down a column), and access
  to the underlying column arrays.
* Data access by row/col offset or label.  Labels are reasonably well
  optimized for both numeric lables and arbitrary (eg, sting) labels.

It does not currently support:
* Views on matrix subset - for currently known access patterns,
  it is more effiicent to copy on subsetting, optimizing for access
  speed over memory use.
* JS iterators - they are too slow.  Use explicit iteration over
  offest or labels.

Important assumptions embedded in the API:
* Columns are implicitly categorical if they are a JS Array and numeric
  (aka continuous) if they are a TypedArray.

There are three index types for row/col indexing:
* IdentityInt32Index - noop index, where the index label is the offset.
* KeyIndex - index arbitrary JS objects.
* DenseInt32Index - integer indexing.  Optimization over KeyIndex as it uses
  Int32Array as a back-map to offsets.  This means that the index array
  must be sized to [minLabel, maxLabel), so this is only useful when the label
  range is relatively close the underlying offset range [minOffset, maxOffset).

All private functions/methods/fields are prefixed by '__', eg, __compile().
Don't use them outside of this file.

Simple example:

  // default indexing is integer offset.
  const df = Dataframe.create([2,2], [['a', 'b'], [0, 1]])
  console.log(df.at(0,0));  // outputs: a
  console.log(df.col(1).asArray());   // outputs: [0, 1]

  // KeyIndex
  const df = new Dataframe([1,2], [['a'], ['b']], null, new KeyIndex(['A', 'B']))
  console.log(df.at(0, 'A')); // outputs: a
  console.log(df.col('A').asArray(); // outputs: ['a']

Performance tuning is primarily focused on columnar access patterns, which is the
dominant pattern in cellxgene.
*/

/**
Dataframe
**/

class Dataframe {
  // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
  __columns: any;

  // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
  __columnsAccessor: any;

  // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
  __id: any;

  // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
  colIndex: any;

  // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
  dims: any;

  // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
  length: any;

  // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
  rowIndex: any;
  /**
  Constructors & factories
  **/

  constructor(
    // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any -- - FIXME: disabled temporarily on migrate to TS.
    dims: any,
    // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any -- - FIXME: disabled temporarily on migrate to TS.
    columnarData: any,
    rowIndex = null,
    colIndex = null,
    __columnsAccessor = [] // private interface
  ) {
    /*
    The base constructor is relatively hard to use - as an alternative,
    see factory methods and clone/slice, below.

    Parameters:
      * dims - 2D array describing intendend dimensionality: [nRows,nCols].
      * columnarData - JS array, nCols in length, containing array
        or TypedArray of length nRows.
      * rowIndex/colIndex - null (create default index using offsets as key),
        or a caller-provided index.
      * __columnsAccessor - private interface, do not specify.  Used internally
        to improve caching of column accessors when possible (eg, clone(), 
        dropCol(), withCol()).
    All columns and indices must have appropriate dimensionality.
    */
    const [nRows, nCols] = dims;
    if (nRows < 0 || nCols < 0) {
      throw new RangeError("Dataframe dimensions must be positive");
    }
    if (!rowIndex) {
      // @ts-expect-error ts-migrate(2322) FIXME: Type 'IdentityInt32Index' is not assignable to typ... Remove this comment to see the full error message
      rowIndex = new IdentityInt32Index(nRows);
    }
    if (!colIndex) {
      // @ts-expect-error ts-migrate(2322) FIXME: Type 'IdentityInt32Index' is not assignable to typ... Remove this comment to see the full error message
      colIndex = new IdentityInt32Index(nCols);
    }
    Dataframe.__errorChecks(dims, columnarData, rowIndex, colIndex);

    this.__columns = Array.from(columnarData);
    this.dims = dims;
    this.length = nRows; // convenience accessor for row dimension
    this.rowIndex = rowIndex;
    this.colIndex = colIndex;
    this.__id = __getMemoId();

    this.__compile(__columnsAccessor);
    Object.freeze(this);
  }

  // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types --- FIXME: disabled temporarily on migrate to TS.
  static __errorChecks(
    // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any -- - FIXME: disabled temporarily on migrate to TS.
    dims: any,
    // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any -- - FIXME: disabled temporarily on migrate to TS.
    columnarData: any,
    // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any -- - FIXME: disabled temporarily on migrate to TS.
    rowIndex: any,
    // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any -- - FIXME: disabled temporarily on migrate to TS.
    colIndex: any
  ) {
    const [nRows, nCols] = dims;

    /* check for expected types */
    if (!Array.isArray(columnarData)) {
      throw new TypeError("Dataframe constructor requires array of columns");
    }
    if (!columnarData.every((c) => isAnyArray(c))) {
      throw new TypeError("Dataframe columns must all be Array or TypedArray");
    }
    if (!isLabelIndex(rowIndex)) {
      throw new TypeError("Dataframe rowIndex is an unsupported type.");
    }
    if (!isLabelIndex(colIndex)) {
      throw new TypeError("Dataframe colIndex is an unsupported type.");
    }

    /* check for expected dimensionality / size */
    if (
      nCols !== columnarData.length ||
      !columnarData.every((c) => c.length === nRows)
    ) {
      throw new RangeError(
        "Dataframe dimension does not match provided data shape"
      );
    }
    if (nRows !== rowIndex.size()) {
      throw new RangeError(
        "Dataframe rowIndex must have same size as underlying data"
      );
    }
    if (nCols !== colIndex.size()) {
      throw new RangeError(
        "Dataframe colIndex must have same size as underlying data"
      );
    }
  }

  // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any -- - FIXME: disabled temporarily on migrate to TS.
  static __compileColumn(column: any, getRowByOffset: any, getRowByLabel: any) {
    /*
      Each column accessor is a function which will lookup data by
      index (ie, is equivalent to dataframe.get(row, col), where 'col'
      is fixed.

      In addition, each column accessor has several functions:

      asArray() -- return the entire column as a native Array or TypedArray.
        Crucially, this native array only supports label indexing.
        Example:
          const arr = df.col('a').asArray();

      has(rlabel) -- return boolean indicating of the row label
        is contained within the column.  Example:
          const isInColumn = df.col('a').includes(99)
        For the default offset indexing, this is identical to:
          const isInColumn = (99 > 0) && (99 < df.nRows);

      ihas(roffset) -- same as has(), but accepts a row offset
        instead of a row label.

      indexOf(value) -- return the label (not offset) of the first instance of
        'value' in the column.  If you want the offset, just use the builtin JS
        indexOf() function, available on both Array and TypedArray.

      iget(offset) -- return the value at 'offset'

      ... and more ...

    */
    const { length } = column;
    const __id = __getMemoId();

    /* get value by row label */
    // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
    const get = function get(rlabel: any) {
      return column[getRowByOffset(rlabel)];
    };

    /* get value by row offset */
    // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
    const iget = function iget(roffset: any) {
      return column[roffset];
    };

    /* full column array access */
    const asArray = function asArray() {
      return column;
    };

    /* test for row label inclusion in column */
    // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
    const has = function has(rlabel: any) {
      const offset = getRowByOffset(rlabel);
      return offset >= 0 && offset < length;
    };

    // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
    const ihas = function ihas(offset: any) {
      return offset >= 0 && offset < length;
    };

    /*
    return first label (index) at which the value is found in this column,
    or undefined if not found.

    NOTE: not found return is DIFFERENT than the default Array.indexOf as
    -1 is a plausible Dataframe row/col label.
    */
    // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
    const indexOf = function indexOf(value: any) {
      const offset = column.indexOf(value);
      if (offset === -1) {
        return undefined;
      }
      return getRowByLabel(offset);
    };

    /*
    Summarize the column data. Lazy eval, memoized
    */
    const summarizeCategorical = callOnceLazy(() =>
      _summarizeCategorical(column)
    );
    const summarize = callOnceLazy(() =>
      isTypedArray(column)
        ? summarizeContinuous(column)
        : summarizeCategorical(column)
    );

    /*
    Create histogram bins for this column.  Memoized.
    */
    const _memoHistoCat = memoize(_histogramCategorical, hashCategorical);
    // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
    const histogramCategorical = (by: any) => _memoHistoCat(get, by);
    let histogram = null;
    if (isTypedArray(column)) {
      const mFn = memoize(histogramContinuous, hashContinuous);
      // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
      histogram = (bins: any, domain: any, by: any) =>
        mFn(get, bins, domain, by);
    } else {
      histogram = histogramCategorical;
    }

    get.summarize = summarize;
    get.summarizeCategorical = summarizeCategorical;
    get.histogram = histogram;
    get.histogramCategorical = histogramCategorical;
    get.asArray = asArray;
    get.has = has;
    get.ihas = ihas;
    get.indexOf = indexOf;
    get.iget = iget;
    get.__id = __id;

    Object.freeze(get);
    return get;
  }

  // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any -- - FIXME: disabled temporarily on migrate to TS.
  __compile(accessors: any) {
    /*
    Compile data accessors for each column.

    Use an existing accessor if provided, else compile a new one.
    */
    const getRowByOffset = this.rowIndex.getOffset.bind(this.rowIndex);
    const getRowByLabel = this.rowIndex.getLabel.bind(this.rowIndex);
    // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
    this.__columnsAccessor = this.__columns.map((column: any, idx: any) => {
      if (accessors[idx]) {
        return accessors[idx];
      }
      return Dataframe.__compileColumn(column, getRowByOffset, getRowByLabel);
    });
    Object.freeze(this.__columnsAccessor);
  }

  // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types --- FIXME: disabled temporarily on migrate to TS.
  clone() {
    /*
    Clone this dataframe
    */
    // @ts-expect-error ts-migrate(2351) FIXME: This expression is not constructable.
    return new this.constructor(
      this.dims,
      [...this.__columns],
      this.rowIndex,
      this.colIndex,
      [...this.__columnsAccessor]
    );
  }

  // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any -- - FIXME: disabled temporarily on migrate to TS.
  withCol(label: any, colData: any, withRowIndex = null) {
    /*
    Create a new DF, which is `this` plus the new column. Example:
    const newDf = df.withCol("foo", [1,2,3]);

    Dimensionality of new column must match existing dataframe.

    Special case: empty dataframe will accept any size column.  Example:
    const newDf = Dataframe.empty().withCol("foo", [1,2,3]);

    If `withRowIndex` specified, the provided index will become the
    rowIndex for the newly created dataframe.   If not specified,
    the rowIndex from `this` will be used (ie, the rowIndex is
    unchanged).
    */
    let dims;
    let rowIndex;
    if (this.isEmpty()) {
      dims = [colData.length, 1];
      rowIndex = null;
    } else {
      dims = [this.dims[0], this.dims[1] + 1];
      ({ rowIndex } = this);
    }

    if (withRowIndex) {
      rowIndex = withRowIndex;
    }

    const columns = [...this.__columns];
    columns.push(colData);
    const colIndex = this.colIndex.withLabel(label);
    const columnsAccessor = [...this.__columnsAccessor];
    // @ts-expect-error ts-migrate(2351) FIXME: This expression is not constructable.
    return new this.constructor(
      dims,
      columns,
      rowIndex,
      colIndex,
      columnsAccessor
    );
  }

  // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any -- - FIXME: disabled temporarily on migrate to TS.
  withColsFrom(dataframe: any, labels: any) {
    /*
    return a new dataframe containing all columns from both `this` and the
    provided dataframe argument.

    The row index from `this` will be used.  All dataframes must have identical
    dimensionality, and no overlapping columns labels.

    Special case, if either dataframe is empty, the other is returned unchanged.

    Arguments:
    * dataframe: a dataframe to combine with `this`
    * labels: columns to pull from `dataframe` and combine with `this`.  If falsey,
      all columns are used.  If an array, must contain a list of labels.  If an
      Object or Map, the key is the columns to pull, which will be stored into the
      new dataframe as the value.

    Example:

    newDf = df.withColsFrom(otherDf);         // combines all columns from both
    newDf = df.withColsFrom(otherDf, ['a']);  // combines df with otherDf['a']
    newDf = df.withColsFrom(otherDf, {a: 'b'}); // combines df with otherDf['a'], but calls it 'b'

     */

    // resolve the source and dest label names.
    let srcLabels;
    let dstLabels;
    if (!labels) {
      // combine all columns
      dstLabels = dataframe.colIndex.labels();
      srcLabels = dstLabels;
    } else if (Array.isArray(labels)) {
      // combine subset of keys with no aliasing
      dstLabels = labels;
      srcLabels = labels;
    } else if (labels instanceof Map) {
      // aliasing with a Map
      srcLabels = Array.from(labels.keys());
      dstLabels = Array.from(labels.values());
    } else {
      // aliasing with an Object
      srcLabels = Object.keys(labels);
      dstLabels = Object.values(labels);
    }

    // if datafame is empty, and no specific labels specified, noop.
    if (dataframe.isEmpty()) {
      if (!labels || srcLabels.length === 0) return this;
      throw new Error("Empty dataframe, unable to pick columns");
    }

    if (this.isEmpty()) {
      // 1. subset dataframe from source keys
      // 2. alias names
      dataframe = dataframe.subset(null, srcLabels);
      for (let i = 0; i < srcLabels.length; i += 1) {
        dataframe = dataframe.renameCol(srcLabels[i], dstLabels[i]);
      }
      return dataframe;
    }

    // otherwise, bulid a new dataframe combining columns from both

    // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
    const srcOffsets = srcLabels.map((l: any) =>
      dataframe.colIndex.getOffset(l)
    );

    // check for label collisions
    if (dstLabels.some(this.hasCol, this)) {
      throw new Error("duplicate key collision");
    }

    // const dims = [this.dims[0], this.dims[1] + dataframe.dims[1]];
    const dims = [this.dims[0], this.dims[1] + srcOffsets.length];
    const { rowIndex } = this;
    const columns = [
      ...this.__columns,
      // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
      ...srcOffsets.map((i: any) => dataframe.__columns[i]),
    ];
    const colIndex = this.colIndex.withLabels(dstLabels);
    const columnsAccessor = [
      ...this.__columnsAccessor,
      // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
      ...srcOffsets.map((i: any) => dataframe.__columnsAccessor[i]),
    ];

    // @ts-expect-error ts-migrate(2351) FIXME: This expression is not constructable.
    return new this.constructor(
      dims,
      columns,
      rowIndex,
      colIndex,
      columnsAccessor
    );
  }

  // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types --- FIXME: disabled temporarily on migrate to TS.
  withColsFromAll(dataframes = []) {
    dataframes = Array.isArray(dataframes) ? dataframes : [dataframes];
    // @ts-expect-error ts-migrate(2554) FIXME: Expected 2 arguments, but got 1.
    return dataframes.reduce((acc, df) => acc.withColsFrom(df), this);
  }

  // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any -- - FIXME: disabled temporarily on migrate to TS.
  dropCol(label: any) {
    /*
    Create a new dataframe, omitting one columns.

    const newDf = df.dropCol("colors");

    Corner case to manage: if dropping the last column, return an empty dataframe.
    */
    if (!this.hasCol(label)) {
      throw new RangeError(`unknown label: ${label}`);
    }

    /* 
    Corner case to manage: if dropping the last column, return an empty dataframe. 
    */
    if (this.dims[1] === 1) {
      return Dataframe.empty();
    }

    const dims = [this.dims[0], this.dims[1] - 1];
    const coffset = this.colIndex.getOffset(label);
    const columns = [...this.__columns];
    columns.splice(coffset, 1);
    const colIndex = this.colIndex.dropLabel(label);
    const columnsAccessor = [...this.__columnsAccessor];
    columnsAccessor.splice(coffset, 1);
    // @ts-expect-error ts-migrate(2351) FIXME: This expression is not constructable.
    return new this.constructor(
      dims,
      columns,
      this.rowIndex,
      colIndex,
      columnsAccessor
    );
  }

  // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any -- - FIXME: disabled temporarily on migrate to TS.
  renameCol(oldLabel: any, newLabel: any) {
    /*
    Accelerator for dropping a column and then adding it again with a new label
    */
    const coffset = this.colIndex.getOffset(oldLabel);
    const colIndex = this.colIndex.dropLabel(oldLabel).withLabel(newLabel);

    const columns = [...this.__columns];
    columns.push(columns[coffset]);
    columns.splice(coffset, 1);

    const columnsAccessor = [...this.__columnsAccessor];
    columnsAccessor.push(columnsAccessor[coffset]);
    columnsAccessor.splice(coffset, 1);

    // @ts-expect-error ts-migrate(2351) FIXME: This expression is not constructable.
    return new this.constructor(
      this.dims,
      columns,
      this.rowIndex,
      colIndex,
      columnsAccessor
    );
  }

  // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any -- - FIXME: disabled temporarily on migrate to TS.
  replaceColData(label: any, newColData: any) {
    /*
    Accelerator for dropping a column then adding it again with same
    label and different values.
    */
    const coffset = this.colIndex.getOffset(label);
    const columns = [...this.__columns];
    columns[coffset] = newColData;
    const columnsAccessor = [...this.__columnsAccessor];
    columnsAccessor[coffset] = null;

    // @ts-expect-error ts-migrate(2351) FIXME: This expression is not constructable.
    return new this.constructor(
      this.dims,
      columns,
      this.rowIndex,
      this.colIndex,
      columnsAccessor
    );
  }

  // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types --- FIXME: disabled temporarily on migrate to TS.
  static empty(rowIndex = null, colIndex = null) {
    const dims = [
      // @ts-expect-error ts-migrate(2531) FIXME: Object is possibly 'null'.
      rowIndex ? rowIndex.size() : 0,
      // @ts-expect-error ts-migrate(2531) FIXME: Object is possibly 'null'.
      colIndex ? colIndex.size() : 0,
    ];
    if (dims[0] && dims[1]) throw new Error("not an empty dataframe");
    return new Dataframe(dims, new Array(dims[1]), rowIndex, colIndex);
  }

  // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any -- - FIXME: disabled temporarily on migrate to TS.
  static create(dims: any, columnarData: any) {
    /*
    Create a dataframe from raw columnar data.  All column arrays
    must have the same length.   Identity indexing will be used.

    Example:
    const df = Dataframe.create([2,2], [new Uint32Array(2), new Float32Array(2)]);
    */
    return new Dataframe(dims, columnarData, null, null);
  }

  // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any -- - FIXME: disabled temporarily on migrate to TS.
  __subset(newRowIndex: any, newColIndex: any) {
    const dims = [...this.dims];

    /* subset columns */
    let { __columns, colIndex, __columnsAccessor } = this;
    if (newColIndex) {
      const colOffsets = this.colIndex.getOffsets(newColIndex.labels());
      __columns = new Array(colOffsets.length);
      __columnsAccessor = new Array(colOffsets.length);
      for (let i = 0, l = colOffsets.length; i < l; i += 1) {
        __columns[i] = this.__columns[colOffsets[i]];
        __columnsAccessor[i] = this.__columnsAccessor[colOffsets[i]];
      }
      colIndex = newColIndex;
      dims[1] = colOffsets.length;
    }

    let { rowIndex } = this;
    if (newRowIndex) {
      const rowOffsets = this.rowIndex.getOffsets(newRowIndex.labels());
      // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
      __columns = __columns.map((col: any) => {
        const newCol = new col.constructor(rowOffsets.length);
        for (let i = 0, l = rowOffsets.length; i < l; i += 1) {
          newCol[i] = col[rowOffsets[i]];
        }
        return newCol;
      });
      rowIndex = newRowIndex;
      dims[0] = rowOffsets.length;
      __columnsAccessor = []; // force a recompile
    }

    if (dims[0] === 0 || dims[1] === 0) return Dataframe.empty();
    return new Dataframe(
      dims,
      __columns,
      rowIndex,
      colIndex,
      __columnsAccessor
    );
  }

  // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any -- - FIXME: disabled temporarily on migrate to TS.
  subset(rowLabels: any, colLabels = null, withRowIndex = null) {
    /*
    Subset by row/col labels.

    withRowIndex allows subset with an index, rather than rowLabels.
    If withRowIndex is specified, rowLabels is ignored.
    */
    let rowIndex = null;
    if (withRowIndex) {
      rowIndex = withRowIndex;
    } else if (rowLabels) {
      rowIndex = this.rowIndex.subset(rowLabels);
    }

    let colIndex = null;
    if (colLabels) {
      colIndex = this.colIndex.subset(colLabels);
    }

    return this.__subset(rowIndex, colIndex);
  }

  // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any -- - FIXME: disabled temporarily on migrate to TS.
  isubset(rowOffsets: any, colOffsets = null, withRowIndex = null) {
    /*
    Subset by row/col offset.

    withRowIndex allows assignment of new row index during subset operation.
    If withRowIndex === null, it will reset the index to identity (offset)
    indexing. If withRowIndex is a label index object, it will be used
    for the new dataframe.
    */
    let rowIndex = null;
    if (withRowIndex) {
      rowIndex = withRowIndex;
    } else if (rowOffsets) {
      rowIndex = this.rowIndex.isubset(rowOffsets);
    }

    let colIndex = null;
    if (colOffsets) {
      colIndex = this.colIndex.isubset(colOffsets);
    }

    return this.__subset(rowIndex, colIndex);
  }

  // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any -- - FIXME: disabled temporarily on migrate to TS.
  isubsetMask(rowMask: any, colMask = null, withRowIndex = null) {
    /*
    Subset on row/column based upon a truthy/falsey array (a mask).

    withRowIndex allows assignment of new row index during subset operation.
    If withRowIndex === null, it will reset the index to identity (offset)
    indexing.  if withRowIndex is a label index object, it will be used
    for the new dataframe.
    */
    const [nRows, nCols] = this.dims;
    if (
      (rowMask && rowMask.length !== nRows) ||
      // @ts-expect-error ts-migrate(2531) FIXME: Object is possibly 'null'.
      (colMask && colMask.length !== nCols)
    ) {
      throw new RangeError("boolean arrays must match row/col dimensions");
    }

    /* convert masks to lists - method wastes space, but is fast */
    // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
    const toList = (mask: any, maxSize: any) => {
      if (!mask) {
        return null;
      }
      const list = new Int32Array(maxSize);
      let elems = 0;
      for (let i = 0, l = mask.length; i < l; i += 1) {
        if (mask[i]) {
          list[elems] = i;
          elems += 1;
        }
      }
      return new Int32Array(list.buffer, 0, elems);
    };
    const rowOffsets = toList(rowMask, nRows);
    const colOffsets = toList(colMask, nCols);
    // @ts-expect-error ts-migrate(2345) FIXME: Argument of type 'Int32Array | null' is not assign... Remove this comment to see the full error message
    return this.isubset(rowOffsets, colOffsets, withRowIndex);
  }

  /**
  Data access with row/col.
  **/

  // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types --- FIXME: disabled temporarily on migrate to TS.
  columns() {
    /* return all column accessors as an array, in offset order */
    return [...this.__columnsAccessor];
  }

  // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any -- - FIXME: disabled temporarily on migrate to TS.
  col(columnLabel: any) {
    /*
    Return accessor bound to a column.  Allows random row access
    based upon the row indexing.  Returns undefined if the
    columnLabel is not present in the dataframe.

    Example for a dataframe with string labeled columns, and
    default (offset) indices for rows (eg, [0, 'foo'])

      const getValue = df.col('foo');
      for (let r = 0; r < df.nRows; r += 1) {
        console.log(r, getValue(r));
      }

    See __compile() for the functions available in a column accessor.
    */
    const coff = this.colIndex.getOffset(columnLabel);
    return this.__columnsAccessor[coff];
  }

  // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any -- - FIXME: disabled temporarily on migrate to TS.
  icol(columnOffset: any) {
    /*
    Return column accessor by offset.
    */
    return Number.isInteger(columnOffset)
      ? this.__columnsAccessor[columnOffset]
      : undefined;
  }

  // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any -- - FIXME: disabled temporarily on migrate to TS.
  at(r: any, c: any) {
    /*
    Access a single value, for a row/col label pair.

    For performance reasons, there are no bounds or existance
    checks on labels, and no defined behavior when these are supplied.
    May return undefined, throw an Error, or do something else for
    non-existant labels.  If you want predictable out-of-bounds
    behavior, use has(), eg,

    const myVal = df.has(r,l) ? df.at(r,l) : undefined;
    */
    const coff = this.colIndex.getOffset(c);
    const roff = this.rowIndex.getOffset(r);
    return this.__columns[coff][roff];
  }

  // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any -- - FIXME: disabled temporarily on migrate to TS.
  iat(r: any, c: any) {
    /*
    Access a single value, for a row/col offset (integer) position.

    For performance reasons, there are no bounds checks on row/col offsets
    or other well-defined behavior for out-of-bounds values.  If you want
    well-defined bounds checking, use ihas(), eg,

    const myVal = df.ihas(r, c) ? df.iat(r, c) : undefined;
    */
    return this.__columns[c][r];
  }

  // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any -- - FIXME: disabled temporarily on migrate to TS.
  has(r: any, c: any) {
    /*
    Test if row/col labels exist in the dataframe - returns true/false
    */
    const [nRows, nCols] = this.dims;
    const coff = this.colIndex.getOffset(c);
    const roff = this.rowIndex.getOffset(r);
    return coff >= 0 && coff < nCols && roff >= 0 && roff < nRows;
  }

  // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any -- - FIXME: disabled temporarily on migrate to TS.
  ihas(r: any, c: any) {
    /*
    Test if row/col offset (integer) position exists in the
    dataframe - returns true/false
    */
    const [nRows, nCols] = this.dims;
    return (
      Number.isInteger(r) &&
      Number.isInteger(c) &&
      c >= 0 &&
      c < nCols &&
      r >= 0 &&
      r < nRows
    );
  }

  // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any -- - FIXME: disabled temporarily on migrate to TS.
  hasCol(c: any) {
    /*
    Test if col label exists - return true/false
    */
    return !!this.col(c);
  }

  // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types --- FIXME: disabled temporarily on migrate to TS.
  isEmpty() {
    /*
    Return true if this is an empty dataframe, ie, has dimensions [0,0]
    */
    const [rows, cols] = this.dims;
    return rows === 0 || cols === 0;
  }

  /****
  Functional (map/reduce/etc) data access

  TODO: most are not yet implemented, as there is no clear use case.   Can easily
  add these as useful.
  ****/

  // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any -- - FIXME: disabled temporarily on migrate to TS.
  mapColumns(callback: any) {
    /*
    map all columns in the dataframe, returning a new dataframe comprised of the
    return values, with the same index as the original dataframe.

    callback MUST not modify the column, but instead return a mutated copy.
    */
    // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
    const columns = this.__columns.map((colData: any, colIdx: any) =>
      callback(colData, colIdx, this)
    );
    // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
    const columnsAccessor = columns.map((c: any, idx: any) =>
      this.__columns[idx] === c ? this.__columnsAccessor[idx] : undefined
    );
    // @ts-expect-error ts-migrate(2351) FIXME: This expression is not constructable.
    return new this.constructor(
      this.dims,
      columns,
      this.rowIndex,
      this.colIndex,
      columnsAccessor
    );
  }

  /*
  Map & reduce of column or row

  TODO remainder of map/reduce functions:  mapCol, mapRow, reduceRow, ...
  */
  /* comment out until we have a use for this

  reduceCol(clabel, callback, initialValue) {
    const coff = this.colIndex.getOffset(clabel);
    const column = this.__columns[coff];
    let start = 0;
    let acc = initialValue;
    if (initialValue === undefined) {
      acc = column[0];
      start = 1;
    }
    for (let i = start, l = column.length; i < l; i += 1) {
      acc = callback(acc, column[i]);
    }
    return acc;
  }
  */
}

export default Dataframe;
