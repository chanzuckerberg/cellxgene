



const joyParser = (data, count = 20) => {
  const genes = [];

  /* setup */

  for (let i = 0; i < count; i++) {
    const gene = {
      key: data.genes[i], /* key values naming: https://bl.ocks.org/armollica/3b5f83836c1de5cca7b1d35409a013e3 */
      values: []
    }

    data.cells.forEach((cell) => {
      gene.values.push({
        value: cell["e"][i]
      })
    })

    genes.push(gene)
  }



  return genes;
}

export default joyParser;
