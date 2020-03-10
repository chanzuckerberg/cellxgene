const maybeTruncateString = (str, maxLength) => {

  let truncatedString = null;
  if (
    str.length > maxLength
  ) {
    truncatedString = `${str.slice(
      0,
      maxLength / 2
    )}…${str.slice(
      -maxLength / 2
    )}`;
  } 

  return truncatedString || str;
}

export default maybeTruncateString;