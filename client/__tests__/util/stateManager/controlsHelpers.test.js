/*
test controls helpers
*/
import { subsetAndResetGeneLists } from "../../../src/util/stateManager/controlsHelpers";
import * as globals from "../../../src/globals";

describe("controls helpers", () => {
  test("subsetAndResetGeneLists", () => {
    const geneList = [];
    const genRandGene = () => Math.random().toString(36).substring(2, 6);

    // build a unique set of genes
    for (let i = 0; i < 150; i += 1) {
      let randGene = genRandGene();
      while (geneList.includes(randGene)) randGene = genRandGene();
      geneList.push(randGene);
    }
    // insert duplicates
    geneList[0] = "dupl";
    geneList[20] = "dupl";
    const state = {
      userDefinedGenes: geneList.slice(0, 20),
      diffexpGenes: geneList.slice(20),
    };
    const [newUserDefinedGenes, newDiffExpGenes] = subsetAndResetGeneLists(
      state
    );
    const expectedNewUserDefinedGenes = [
      ...geneList.slice(0, 20),
      ...geneList.slice(21),
    ].slice(0, globals.maxGenes);
    expect(globals.maxUserDefinedGenes).toBeLessThan(globals.maxGenes);
    expect(geneList.length).toBeGreaterThan(globals.maxGenes);
    expect(newUserDefinedGenes).toHaveLength(globals.maxGenes);
    expect(newUserDefinedGenes).toStrictEqual(expectedNewUserDefinedGenes);
    expect(newDiffExpGenes).toStrictEqual([]);
  });
});
