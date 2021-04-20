/**
 * The function will split string by comma or space, return unique non-empty strings
 * @param geneString - a string of comma delimited genes
 * @returns an array
 */
import uniq from "lodash.uniq";

export default function parseBulkGeneString(geneString) {
  return _.pull(uniq(geneString.split(/[ ,]+/)), "");
}
