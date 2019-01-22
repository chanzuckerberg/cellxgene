import { sort, sortIndex } from "../typedCrossfilter/sort";

/*
Dataframe is an immutable 2D matrix similiar to Pandas Dataframe,
but (currently) without all of the mathematical support functions.
Data is stored in column-major layout, and each column is unimorphic.

It supports:
* Relatively effiicent creation, cloning and subsetting ("cut")
* Very efficient columnar access (eg, sum down a column), and access
  to the underlying column arrays
* Data access by offset or label

It does not currently support:
* Views on matrix subset - for currently envisioned access patterns,
  it is more effiicent to copy on cut/slice
* JS iterators - they are too slow.  Use explicit iteration over
  offest or labels.

All private functions/methods are prefixed by '__', eg, __compile().
*/

/**
Label indexing
**/

/* eslint-disable class-methods-use-this */
class IdentityIndex {
  getOffset(i) {
    return i;
  }

  cut(labelArray) {
    return new IntIndex(labelArray);
  }
}

class RangeIndex {
  constructor(start, stop = 0, step = 1) {
    this.start = start;
    this.stop = stop;
    this.step = step;
    this.__compile();
  }

  __compile() {
    const [start, stop, step] = this;
    const max = Math.trunc((stop - start) / step);
    /* eslint-disable-next-line no-new-func */
    this.getOffset = new Function(
      "i",
      `
        const loc = ${start} + i*${step};
        if (loc > ${max}) {
          throw new RangeError();
        }
        return loc;
      `
    );
  }

  cut(labelArray) {
    return new IntIndex(labelArray);
  }
}

class IntIndex {
  constructor(index) {
    this.index = index;
    this.__compile();
  }

  __compile() {
    this.getOffset = (function getOffsetClosure(index) {
      return function getOffset(i) {
        return index[i];
      };
    })(this.index);
  }

  cut(labelArray) {
    return new IntIndex(labelArray);
  }
}

class KeyIndex {
  constructor(labels) {
    const index = new Map();
    labels.forEach((v, i) => index.set(v, i));
    this.index = index;
    this.__compile();
  }

  __compile() {
    const { index } = this;
    this.getOffset = function getOffset(k) {
      return index.get(k);
    };
  }

  cut(labelArray) {
    return new KeyIndex(labelArray);
  }
}
/* eslint-enable class-methods-use-this */

/**
Dataframe
**/

class Dataframe {
  /**
  Constructors & factories
  **/

  constructor(dims, columnarData, rowLabelIndex, colLabelIndex) {
    /*
    The base constructor is relatively hard to use - as an alternative,
    see factory methods and clone/slice, below.
    */
    Dataframe.__errorChecks(dims, columnarData);
    const [nRows, nCols] = dims;
    if (!rowLabelIndex) {
      rowLabelIndex = new IdentityIndex();
    }
    if (!colLabelIndex) {
      colLabelIndex = new IdentityIndex();
    }

    this.columns = Array.from(columnarData);
    this.nRows = nRows;
    this.nCols = nCols;
    this.rowLabelIndex = rowLabelIndex;
    this.colLabelIndex = colLabelIndex;

    this.__compile();
  }

  static __errorChecks(dims, columnarData) {
    // TODO: more error checking on inputs, eg,
    //    - column type matches schema, i.e. array of T
    //    - column names are unique, and legal variable names
    //    - fewer than max (??) columns

    const [nRows, nCols] = dims;
    if (nRows <= 0 || nCols <= 0) {
      throw new RangeError("Dataframe dimensions must be positive");
    }
    if (
      nCols !== columnarData.length ||
      !columnarData.every(c => c.length === nRows)
    ) {
      throw new RangeError(
        "Dataframe constructor parameters do not have the same shape"
      );
    }
  }

  __compile() {
    /*
    Compile data accessors for each column.
    */
    const roff = this.rowLabelIndex.getOffset;
    this.columnsAccessor = this.columns.map(column => {
      const { length } = column;

      /* default value access by row label */
      const getCol = function getCol(rlabel) {
        return column[roff(rlabel)];
      };

      /* full column array access */
      const asArray = function asArray() {
        return column;
      };

      /* test for row label inclusion in column */
      const includes = function includes(rlabel) {
        const offset = roff(rlabel);
        return offset >= 0 && offset < length;
      };

      getCol.asArray = asArray;
      getCol.includes = includes;
      return getCol;
    });
  }

  clone() {
    return new this.constructor(
      [this.nRows, this.nCols],
      [...this.columns],
      this.rowLabelIndex,
      this.colLabelIndex
    );
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

  cut(rowLabels, colLabels) {
    /*
    Cut (aka slice) a dataframe by row/column label.  Accepted values are:
      - array of labels indicating the desired row and/or columns
      - falsey value, indicating that this dimension should not be cut
    Example:
      const newDf = df.cut([0,2,4], null);
    */

    /* algo:

    default label arrays (may be falsey/all)
    convert label arrays to offsets (may be falsey/all)
    set dims
    cut column & columnLabelIndex
    cut rows & rowLabelIndex
    call constructor

    */
    const dims = [this.nRows, this.nCols];

    const getSortedLabelAndOffsets = function getSortedLabelAndOffsets(
      labels,
      index
    ) {
      if (!labels || !index) {
        return [null, null];
      }
      const offsets = labels.map(label => {
        const off = index.getOffset(label);
        if (off === undefined) {
          throw new RangeError(`unknown label: ${label}`);
        }
        return off;
      });
      const sortedLabels = sortIndex(labels, offsets);
      const sortedOffsets = sort(offsets);
      return [sortedLabels, sortedOffsets];
    };

    let colOffsets;
    let { colLabelIndex } = this;
    if (colLabels) {
      [colLabels, colOffsets] = getSortedLabelAndOffsets(
        colLabels,
        this.colLabelIndex
      );
      dims[1] = colLabels.length;
      colLabelIndex = this.colLabelIndex.cut(colLabels);
    }

    let rowOffsets;
    let { rowLabelIndex } = this;
    if (rowLabels) {
      [rowLabels, rowOffsets] = getSortedLabelAndOffsets(
        rowLabels,
        this.rowLabelIndex
      );
      dims[0] = rowLabels.length;
      rowLabelIndex = this.rowLabelIndex.cut(rowLabels);
    }

    /* cut columns */
    let { columns } = this;
    if (colLabels) {
      columns = colOffsets.map(coff => this.columns[coff]);
    }

    /* cut rows */
    if (rowLabels) {
      columns = columns.map(col => {
        const newCol = new col.constructor(rowOffsets.length);
        for (let i = 0, l = rowOffsets.length; i < l; i += 1) {
          newCol[i] = col[rowOffsets[i]];
        }
        return newCol;
      });
    }

    return new Dataframe(dims, columns, rowLabelIndex, colLabelIndex);
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

    In addition, each bound column accessor has several more functions:

    asArray() -- return the entire column as a native Array or TypedArray.
      Crucially, this native array only supports label indexing.
      Example:
        const arr = df.col('a').asArray();

    includes(rlabel) -- return boolean indicating of the row label
      is contained within the column.  Example:
        const isInColumn = df.col('a').includes(99)
      For the default offset indexing, this is identical to:
        const isInColumn = (99 > 0) && (99 < df.nRows);
    */
    const coff = this.colLabelIndex.getOffset(columnLabel);
    return this.columnsAccessor[coff];
  }

  row(rowLabel) {
    /*
    Return accessor bound to a row.  Allows random column access.
    Example:

      console.log("Value at [99,'foo']:", df.row(99)("foo"));

    */
    const roff = this.rowLabelIndex.getOffset(rowLabel);
    if (roff === undefined) {
      return undefined;
    }
    const coff = this.colLabelIndex.getOffset;
    const { columns } = this;
    return function col(clabel) {
      return columns[coff(clabel)][roff];
    };
  }

  at(r, c) {
    /*
      Access a single value, for a row/col label pair
    */
    const coff = this.colLabelIndex.getOffset(c);
    const roff = this.rowLabelIndex.getOffset(r);
    return this.columns[coff][roff];
  }

  iat(r, c) {
    /*
    Access a single value, for a row/col offset (integer) position
    */
    return this.columns[c][r];
  }

  /****
  Functional (map/reduce/etc) data access
  ****/

  /*
  Map & reduce of column or row

  XXX TODO remainder of map/reduce functions:  mapCol, mapRow, reduceRow, ...
  */
  reduceCol(clabel, callback, initialValue) {
    const coff = this.colLabelIndex.getOffset(clabel);
    const column = this.columns[coff];
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
}

export { Dataframe, RangeIndex, IntIndex, IdentityIndex, KeyIndex };
