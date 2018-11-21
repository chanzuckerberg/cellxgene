// jshint esversion: 6

// values is [ [optVal, optIdx], ...]
// index is range array
// return sorted index
export default values =>
  values.sort((a, b) => {
    const textA = String(a[0]).toUpperCase();
    const textB = String(b[0]).toUpperCase();
    return textA < textB ? -1 : textA > textB ? 1 : 0;
  });
