import { IdentityInt32Index } from "./labelIndex";
// weird cross-dependency that we should clean up someday...
import { sort } from "../typedCrossfilter/sort";
import { isTypedArray, isArrayOrTypedArray, callOnceLazy } from "./util";
import { summarizeContinuous, summarizeCategorical } from "./summarize";

/*
Dataframe is an immutable 2D matrix similiar to Python Pandas Dataframe,
but (currently) without all of the surrounding support functions.
Data is stored in column-major layout, and each column is monomorphic.

It supports:
* Relatively efficient creation, cloning and subsetting ("cut")
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
  /**
  Constructors & factories
  **/

  constructor(dims, columnarData, rowIndex = null, colIndex = null) {
    /*
    The base constructor is relatively hard to use - as an alternative,
    see factory methods and clone/slice, below.

    Parameters:
      * dims - 2D array describing intendend dimensionality: [nRows,nCols].
      * columnarData - JS array, nCols in length, containing array
        or TypedArray of length nRows.
      * rowIndex/colIndex - null (create default index using offsets as key),
        or a caller-provided index.
    All columns and indices must have appropriate dimensionality.
    */
    Dataframe.__errorChecks(dims, columnarData, rowIndex, colIndex);
    const [nRows, nCols] = dims;
    if (!rowIndex) {
      rowIndex = new IdentityInt32Index(nRows);
    }
    if (!colIndex) {
      colIndex = new IdentityInt32Index(nCols);
    }

    this.__columns = Array.from(columnarData);
    this.dims = dims;
    this.length = nRows; // convenience accessor for row dimension
    this.rowIndex = rowIndex;
    this.colIndex = colIndex;

    this.__compile();
  }

  static __errorChecks(dims, columnarData) {
    const [nRows, nCols] = dims;
    if (nRows < 0 || nCols < 0) {
      throw new RangeError("Dataframe dimensions must be positive");
    }
    if (!Array.isArray(columnarData)) {
      throw new TypeError("Dataframe constructor requires array of columns");
    }
    if (!columnarData.every(c => isArrayOrTypedArray(c))) {
      throw new TypeError("Dataframe columns must all be Array or TypedArray");
    }
    if (
      nCols !== columnarData.length ||
      !columnarData.every(c => c.length === nRows)
    ) {
      throw new RangeError(
        "Dataframe dimension does not match column data shape"
      );
    }
  }

  __compile() {
    /*
    Compile data accessors for each column.

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

    */
    const { getOffset, getLabel } = this.rowIndex;
    this.__columnsAccessor = this.__columns.map(column => {
      const { length } = column;

      /* get value by row label */
      const get = function get(rlabel) {
        return column[getOffset(rlabel)];
      };

      /* get value by row offset */
      const iget = function iget(roffset) {
        return column[roffset];
      };

      /* full column array access */
      const asArray = function asArray() {
        return column;
      };

      /* test for row label inclusion in column */
      const has = function has(rlabel) {
        const offset = getOffset(rlabel);
        return offset >= 0 && offset < length;
      };

      const ihas = function ihas(offset) {
        return offset >= 0 && offset < length;
      };

      /*
      return first label (index) at which the value is found in this column,
      or undefined if not found.

      NOTE: not found return is DIFFERENT than the default Array.indexOf as
      -1 is a plausible Dataframe row/col label.
      */
      const indexOf = function indexOf(value) {
        const offset = column.indexOf(value);
        if (offset === -1) {
          return undefined;
        }
        return getLabel(offset);
      };

      /*
      Summarize the column data. Lazy eval;
      */
      const summarize = callOnceLazy(() =>
        isTypedArray(column)
          ? summarizeContinuous(column)
          : summarizeCategorical(column)
      );

      get.summarize = summarize;
      get.asArray = asArray;
      get.has = has;
      get.ihas = ihas;
      get.indexOf = indexOf;
      get.iget = iget;
      return get;
    });
  }

  clone() {
    /*
    Clone this dataframe
    */
    return new this.constructor(
      this.dims,
      [...this.__columns],
      this.rowIndex,
      this.colIndex
    );
  }

  withCol(label, colData) {
    /*
    Like clone(), but adds a column.  Example:
    const newDf = df.withCol("foo", [1,2,3]);

    Dimensionality of new column must match existing dataframe.

    Special case: empty dataframe will accept any size column.  Example:
    const newDf = Dataframe.empty().withCol("foo", [1,2,3]);
    */
    let dims;
    if (this.isEmpty()) {
      dims = [colData.length, 1];
    } else {
      dims = [this.dims[0], this.dims[1] + 1];
    }

    const columns = [...this.__columns];
    columns.push(colData);
    const colIndex = this.colIndex.with(label);
    return new this.constructor(dims, columns, this.rowIndex, colIndex);
  }

  dropCol(label) {
    /*
    Create a new dataframe, omitting one columns.

    const newDf = df.dropCol("colors");
    */
    const dims = [this.dims[0], this.dims[1] - 1];
    const coffset = this.colIndex.getOffset(label);
    const columns = [...this.__columns];
    columns.splice(coffset, 1);
    const colIndex = this.colIndex.drop(label);
    return new this.constructor(dims, columns, this.rowIndex, colIndex);
  }

  static empty(rowIndex = null, colIndex = null) {
    return new Dataframe([0, 0], [], rowIndex, colIndex);
  }

  static create(dims, columnarData) {
    /*
    Create a dataframe from raw columnar data.  All column arrays
    must have the same length.   Identity indexing will be used.

    Example:
    const df = Dataframe.create([2,2], [new Uint32Array(2), new Float32Array(2)]);
    */
    return new Dataframe(dims, columnarData, null, null);
  }

  __cut(rowOffsets, colOffsets) {
    const dims = [...this.dims];

    const getSortedLabelAndOffsets = (offsets, index) => {
      /*
      Given offsets, return both offsets and associated lables,
      sorted by offset.
      */
      if (!offsets) {
        return [null, null];
      }
      const sortedOffsets = sort(offsets);
      const sortedLabels = new Array(sortedOffsets.length);
      for (let i = 0, l = sortedOffsets.length; i < l; i += 1) {
        sortedLabels[i] = index.getLabel(sortedOffsets[i]);
      }
      return [sortedLabels, sortedOffsets];
    };

    let { colIndex } = this;
    if (colOffsets) {
      let colLabels;
      [colLabels, colOffsets] = getSortedLabelAndOffsets(
        colOffsets,
        this.colIndex
      );
      dims[1] = colOffsets.length;
      colIndex = this.colIndex.cut(colLabels);
    }

    let { rowIndex } = this;
    if (rowOffsets) {
      let rowLabels;
      [rowLabels, rowOffsets] = getSortedLabelAndOffsets(
        rowOffsets,
        this.rowIndex
      );
      dims[0] = rowLabels.length;
      rowIndex = this.rowIndex.cut(rowLabels);
    }

    /* cut columns */
    let columns = this.__columns;
    if (colOffsets) {
      columns = new Array(colOffsets.length);
      for (let i = 0, l = colOffsets.length; i < l; i += 1) {
        columns[i] = this.__columns[colOffsets[i]];
      }
    }

    /* cut rows */
    if (rowOffsets) {
      columns = columns.map(col => {
        const newCol = new col.constructor(rowOffsets.length);
        for (let i = 0, l = rowOffsets.length; i < l; i += 1) {
          newCol[i] = col[rowOffsets[i]];
        }
        return newCol;
      });
    }
    return new Dataframe(dims, columns, rowIndex, colIndex);
  }

  cutByList(rowLabels, colLabels = null) {
    const toOffsets = (labels, index) => {
      if (!labels) {
        return null;
      }
      return labels.map(label => {
        const off = index.getOffset(label);
        if (off === undefined) {
          throw new RangeError(`unknown label: ${label}`);
        }
        return off;
      });
    };

    const rowOffsets = toOffsets(rowLabels, this.rowIndex);
    const colOffsets = toOffsets(colLabels, this.colIndex);
    return this.__cut(rowOffsets, colOffsets);
  }

  icutByList(rowOffsets, colOffsets = null) {
    return this.__cut(rowOffsets, colOffsets);
  }

  icutByMask(rowMask, colMask = null) {
    /*
    Cut on row/column based upon a truthy/falsey array.
    */
    const [nRows, nCols] = this.dims;
    if (
      (rowMask && rowMask.length !== nRows) ||
      (colMask && colMask.length !== nCols)
    ) {
      throw new RangeError("boolean arrays must match row/col dimensions");
    }

    /* convert masks to lists - method wastes space, but is fast */
    const toList = (mask, maxSize) => {
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
    return this.__cut(rowOffsets, colOffsets);
  }

  /**
  Data access with row/col.
  **/

  col(columnLabel) {
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

  icol(columnOffset) {
    /*
    Return column accessor by offset.
    */
    return this.__columnsAccessor[columnOffset];
  }

  at(r, c) {
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

  iat(r, c) {
    /*
    Access a single value, for a row/col offset (integer) position.

    For performance reasons, there are no bounds checks on row/col offsets
    or other well-defined behavior for out-of-bounds values.  If you want
    well-defined bounds checking, use ihas(), eg,

    const myVal = df.ihas(r, c) ? df.iat(r, c) : undefined;
    */
    return this.__columns[c][r];
  }

  has(r, c) {
    /*
    Test if row/col labels exist in the dataframe - returns true/false
    */
    const [nRows, nCols] = this.dims;
    const coff = this.colIndex.getOffset(c);
    const roff = this.rowIndex.getOffset(r);
    return coff >= 0 && coff < nCols && roff >= 0 && roff < nRows;
  }

  ihas(r, c) {
    /*
    Test if row/col offset (integer) position exists in the
    dataframe - returns true/false
    */
    const [nRows, nCols] = this.dims;
    return c >= 0 && c < nCols && r >= 0 && r < nRows;
  }

  hasCol(c) {
    /*
    Test if col label exists - return true/false
    */
    return !!this.col(c);
  }

  isEmpty() {
    const [rows, cols] = this.dims;
    return rows === 0 && cols === 0;
  }

  /****
  Functional (map/reduce/etc) data access

  XXX: not yet implemented, as there is no clear use case.   Can easily
  add these as useful.
  ****/

  /*
  Map & reduce of column or row

  XXX TODO remainder of map/reduce functions:  mapCol, mapRow, reduceRow, ...
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
