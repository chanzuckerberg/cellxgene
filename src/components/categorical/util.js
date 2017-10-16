export const alphabeticallySortedValues = (values) => {
  return Object.keys(values).sort((a, b) => {
    var textA = a.toUpperCase();
    var textB = b.toUpperCase();
    return (textA < textB) ? -1 : (textA > textB) ? 1 : 0;
  })
}
