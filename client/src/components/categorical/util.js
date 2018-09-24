// jshint esversion: 6
export const alphabeticallySortedValues = values =>
  Object.keys(values).sort((a, b) => {
    const textA = a.toUpperCase();
    const textB = b.toUpperCase();
    return textA < textB ? -1 : textA > textB ? 1 : 0;
  });
