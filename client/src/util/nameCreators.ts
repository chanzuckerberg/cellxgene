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

const makeDimensionName = (namespace: string, key: string): string =>
  `${namespace}_${key}`;

export const layoutDimensionName = (key: string): string =>
  makeDimensionName("layout", key);

export const obsAnnoDimensionName = (key: string): string =>
  makeDimensionName("obsAnno", key);

export const diffexpDimensionName = (key: string): string =>
  makeDimensionName("varData_diffexp", key);

export const userDefinedDimensionName = (key: string): string =>
  makeDimensionName("varData_userDefined", key);

export const geneSetSummaryDimensionName = (key: string): string =>
  makeDimensionName("geneSetSummary", key);

export interface ContinuousNamespace {
  isObs?: boolean;
  isDiffExp?: boolean;
  isUserDefined?: boolean;
  isGeneSetSummary?: boolean;
}

/*

  ie., makeContinuousDimensionName(continuousNamespace = {isObs: true}, "total_reads")
  see: histogram brush, as it doesn't know what type of continuous it was with only field
*/
export const makeContinuousDimensionName = (
  continuousNamespace: ContinuousNamespace,
  key: string
): string => {
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
