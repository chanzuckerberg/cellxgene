/*

The original motivation for this file is name collision. It's possible
to have two of the same gene, user defined and diffexpressed, on the client
simultaneously, looking at the same data in the cache but with two
dimensions and two crossfilters. It's also possible that someone could
have a obsAnnotation named X, but we are using that for layout. So
we namespace, and abstract to avoid proliferating strings throughout the
codebase.

*/

const makeDimensionName = (namespace, key) => `${namespace}_${key}`;

export const layoutDimensionName = (key) => makeDimensionName("layout", key);
export const obsAnnoDimensionName = (key) => makeDimensionName("obsAnno", key);
export const diffexpDimensionName = (key) =>
  makeDimensionName("varData_diffexp", key);
export const userDefinedDimensionName = (key) =>
  makeDimensionName("varData_userDefined", key);

/*
    continuousNamespace = {
    isObs: true,
    isDiffExp: false,
    isUserDefined: false
  }

  ie., makeContinuousDimensionName(continuousNamespace = {isObs: true}, "total_reads")
  see: histogram brush, as it doesn't know what type of continuous it was with only field
*/
export const makeContinuousDimensionName = (continuousNamespace, key) => {
  let name;
  if (continuousNamespace.isObs) {
    name = obsAnnoDimensionName(key);
  } else if (continuousNamespace.isDiffExp) {
    name = diffexpDimensionName(key);
  } else if (continuousNamespace.isUserDefined) {
    name = userDefinedDimensionName(key);
  } else {
    throw new Error("unknown continuous dimension");
  }

  return name;
};
