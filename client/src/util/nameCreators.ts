/*

The original motivation for this file is name collision. It's possible
to have two of the same gene, user defined and diffexpressed, on the client
simultaneously, looking at the same data in the cache but with two
dimensions and two crossfilters. It's also possible that someone could
have a obsAnnotation named X, but we are using that for layout. So
we namespace, and abstract to avoid proliferating strings throughout the
codebase.

It is _no longer_ used to remove collisions in the crossfilter or
anno matrix namespaces. It is still used by the component tier.

*/

// eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
const makeDimensionName = (namespace: any, key: any) => `${namespace}_${key}`;

// eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any -- - FIXME: disabled temporarily on migrate to TS.
export const layoutDimensionName = (key: any) =>
  makeDimensionName("layout", key);
// eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any -- - FIXME: disabled temporarily on migrate to TS.
export const obsAnnoDimensionName = (key: any) =>
  makeDimensionName("obsAnno", key);
// eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any -- - FIXME: disabled temporarily on migrate to TS.
export const diffexpDimensionName = (key: any) =>
  makeDimensionName("varData_diffexp", key);
// eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any -- - FIXME: disabled temporarily on migrate to TS.
export const userDefinedDimensionName = (key: any) =>
  makeDimensionName("varData_userDefined", key);
// eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any -- - FIXME: disabled temporarily on migrate to TS.
export const geneSetSummaryDimensionName = (key: any) =>
  makeDimensionName("geneSetSummary", key);

/*
    continuousNamespace = {
    isObs: true,
    isDiffExp: false,
    isUserDefined: false,
    isGeneSet: false,
  }

  ie., makeContinuousDimensionName(continuousNamespace = {isObs: true}, "total_reads")
  see: histogram brush, as it doesn't know what type of continuous it was with only field
*/
// eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types --- FIXME: disabled temporarily on migrate to TS.
export const makeContinuousDimensionName = (
  // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any -- - FIXME: disabled temporarily on migrate to TS.
  continuousNamespace: any,
  // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any -- - FIXME: disabled temporarily on migrate to TS.
  key: any
) => {
  let name;
  if (continuousNamespace.isObs) {
    name = obsAnnoDimensionName(key);
  } else if (continuousNamespace.isDiffExp) {
    name = diffexpDimensionName(key);
  } else if (continuousNamespace.isUserDefined) {
    name = userDefinedDimensionName(key);
  } else if (continuousNamespace.isGeneSetSummary) {
    name = geneSetSummaryDimensionName(key);
  } else {
    throw new Error("unknown continuous dimension");
  }

  return name;
};
