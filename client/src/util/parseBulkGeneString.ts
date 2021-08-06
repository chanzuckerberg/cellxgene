/**
 * The function will split string by comma or space, return unique non-empty strings
 * @param geneString - a string of comma delimited genes
 * @returns an array
 */
import pull from "lodash.pull";
import uniq from "lodash.uniq";

export default function parseBulkGeneString(geneString: string) {
  return pull(uniq(geneString.split(/[ ,]+/)), "");
}
