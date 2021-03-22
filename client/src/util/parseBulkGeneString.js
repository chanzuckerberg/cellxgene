/**
 * The function will split string by comma or space, return unique non-empty strings
 * @param geneString - a string of comma delimited genes
 * @returns an array
 */
import _ from "lodash";

export default function parseBulkGeneString(geneString) {
  return _.pull(_.uniq(geneString.split(/[ ,]+/)), "");
}
