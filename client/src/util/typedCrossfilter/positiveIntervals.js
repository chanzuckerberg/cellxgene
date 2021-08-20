// Interval operations - very simple version of interval set relationship
// operators.   An interval is a multi-interval list of [min, max),
// where min and max are mandatory.   Constraints:
//    * min <= max, min >= 0
//    * empty interval groups are OK, ie, []
//    * Legal intervals:   [], [ [0, 1], ... ]
//    * Not legal:   [ [] ]
//
// All intervals are represented by simple JS arrays/numbers.
//
// Code assumes intervals have a low cardinality; many operations are done
// with a brute force scan.   Little attempt to reduce GC pressure.
//
class PositiveIntervals {
  // Canonicalize - ensure that:
  //  1. no overlapping intervals
  //  2. sorted in order of interval min.
  //
  static canonicalize(A) {
    if (A.length <= 1) return A;
    const copy = A.slice();
    copy.sort((a, b) => a[0] - b[0]);
    const res = [];
    res.push(copy[0]);
    for (let i = 1, len = copy.length; i < len; i += 1) {
      if (copy[i][0] > res[res.length - 1][1]) {
        // non-overlapping, add to result
        res.push(copy[i]);
      } else if (copy[i][1] > res[res.length - 1][1]) {
        // merge this into previous
        // eslint-disable-line prefer-destructuring -- destructuring impedes readability
        res[res.length - 1][1] = copy[i][1];
      }
    }
    return res;
  }

  // Return interval with values belonging to both A and B. Essentially
  // a set union operation.
  //
  static union(A, B) {
    return PositiveIntervals.canonicalize([...A, ...B]);
  }

  static _flatten(A, B) {
    const points = []; /* point, A, start */
    for (let a = 0; a < A.length; a += 1) {
      points.push([A[a][0], true, true]);
      points.push([A[a][1], true, false]);
    }
    for (let b = 0; b < B.length; b += 1) {
      points.push([B[b][0], false, true]);
      points.push([B[b][1], false, false]);
    }
    // Sort order: point, then start
    points.sort((a, b) => (a[0] !== b[0] ? a[0] - b[0] : a[2] ? 1 : -1));
    return points;
  }

  // A - B, ie, the interval with all values in A that are not in B.  Essentially
  // a set difference operation.
  //
  static difference(A, B) {
    // Corner cases
    if (A.length === 0 || B.length === 0) {
      return PositiveIntervals.canonicalize(A);
    }

    const cA = PositiveIntervals.canonicalize(A);
    const cB = PositiveIntervals.canonicalize(B);

    const points = PositiveIntervals._flatten(cA, cB);
    const res = [];
    let aDepth = 0;
    let depth = 0;
    let intervalStart;
    for (let i = 0; i < points.length; i += 1) {
      const p = points[i];
      const delta = p[2] ? 1 : -1;
      depth += delta;
      if (p[1]) aDepth += delta;

      if (i === points.length - 1 || p[0] !== points[i + 1][0]) {
        if (aDepth === 1 && depth === 1) {
          [intervalStart] = p;
        } else if (intervalStart !== undefined) {
          res.push([intervalStart, p[0]]);
          intervalStart = undefined;
        }
      }
    }
    // guaranteed to be in canonical form
    return res;
  }

  // Return interval with values belonging to A or B.  Essentially a set
  // intersection.
  //
  static intersection(A, B) {
    if (A.length === 0 || B.length === 0) {
      return [];
    }

    const cA = PositiveIntervals.canonicalize(A);
    const cB = PositiveIntervals.canonicalize(B);

    const points = PositiveIntervals._flatten(cA, cB);
    const res = [];
    let depth = 0;
    let intervalStart;
    for (let i = 0; i < points.length; i += 1) {
      const p = points[i];
      depth += p[2] ? 1 : -1;
      if (depth === 2) {
        [intervalStart] = p;
      } else if (intervalStart !== undefined) {
        res.push([intervalStart, p[0]]);
        intervalStart = undefined;
      }
    }
    // guaranteed to be in canonical form
    return res;
  }
}

export default PositiveIntervals;
