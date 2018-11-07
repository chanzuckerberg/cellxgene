const SmallArray = 32;

function insertionsort(a, lo, hi) {
  for (let i = lo + 1; i < hi + 1; i += 1) {
    const x = a[i];
    let j;
    for (j = i; j > lo && a[j - 1] > x; j -= 1) {
      a[j] = a[j - 1];
    }
    a[j] = x;
  }
  return a;
}

function insertionsortIndirect(a, s, lo, hi) {
  for (let i = lo + 1; i < hi + 1; i += 1) {
    const x = a[i];
    const t = s[x];
    let j;
    for (j = i; j > lo && s[a[j - 1]] > t; j -= 1) {
      a[j] = a[j - 1];
    }
    a[j] = x;
  }
  return a;
}

function quicksort(a, lo, hi) {
  if (hi - lo < SmallArray) {
    return insertionsort(a, lo, hi);
  }
  if (lo < hi) {
    // partition
    const mid = Math.floor((lo + hi) / 2);
    const p = a[mid];
    let i = lo - 1;
    let j = hi + 1;
    while (i < j) {
      do {
        i += 1;
      } while (a[i] < p);
      do {
        j -= 1;
      } while (a[j] > p);
      if (i < j) {
        const tmp = a[i];
        a[i] = a[j];
        a[j] = tmp;
      }
    }
    // sort
    quicksort(a, lo, j);
    quicksort(a, j + 1, hi);
  }
  return a;
}

function quicksortIndirect(a, s, lo, hi) {
  if (hi - lo < SmallArray) {
    return insertionsortIndirect(a, s, lo, hi);
  }
  if (lo < hi) {
    // partition
    const mid = Math.floor((lo + hi) / 2);
    const p = a[mid];
    const t = s[p];
    let i = lo - 1;
    let j = hi + 1;
    while (i < j) {
      do {
        i += 1;
      } while (s[a[i]] < t);
      do {
        j -= 1;
      } while (s[a[j]] > t);
      if (i < j) {
        const tmp = a[i];
        a[i] = a[j];
        a[j] = tmp;
      }
    }
    // sort
    quicksortIndirect(a, s, lo, j);
    quicksortIndirect(a, s, j + 1, hi);
  }
  return a;
}

// Convenience wrappers
export function sort(arr, comparator = undefined) {
  if (comparator !== undefined) {
    // XXX for now
    return arr.sort(arr, comparator);
  }
  return quicksort(arr, 0, arr.length - 1);
}

export function sortIndex(index, source) {
  return quicksortIndirect(index, source, 0, index.length - 1);
}
