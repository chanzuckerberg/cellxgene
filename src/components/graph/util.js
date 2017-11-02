createExpressionsCountsMap () {

  const CHANGE_ME_MAGIC_GENE_INDEX = 5;

  const expressionsCountsMap = {};

  /* currently selected gene */
  expressionsCountsMap.geneName = this.state.expressions.data.genes[3];

  let maxExpressionValue = 0;

  /* create map of expressions for every cell */
  this.state.expressions.data.cells.map((c) => {
    /* cellname = 234 */
    expressionsCountsMap[c.cellname] = c["e"][CHANGE_ME_MAGIC_GENE_INDEX];
    /* collect the maximum value as we iterate */
    if (c["e"][CHANGE_ME_MAGIC_GENE_INDEX] > maxExpressionValue) {
      maxExpressionValue = c["e"][CHANGE_ME_MAGIC_GENE_INDEX]
    }
  })

  expressionsCountsMap.maxValue = maxExpressionValue;

  return expressionsCountsMap;
}
