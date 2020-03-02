/*
test controls helpers
*/
import { subsetAndResetGeneLists } from "../../../src/util/stateManager/controlsHelpers";
import * as globals from "../../../src/globals";

describe("controls helpers", () => {

  test("subsetAndResetGeneLists", () => {
    const geneList = [...Array(150).keys()].map(() =>
      Math.random().toString(36).substring(2, 6) // random string of 4 characters
    );
    const state = {
      userDefinedGenes: geneList.slice(0, 20),
      diffexpGenes: geneList.slice(20),
    };
    const [newUserDefinedGenes, newDiffExpGenes] = subsetAndResetGeneLists(state);
    expect(globals.maxUserDefinedGenes).toBeLessThan(globals.maxGenes);
    expect(geneList.length).toBeGreaterThan(globals.maxGenes);
    expect(newUserDefinedGenes).toHaveLength(globals.maxGenes);
    expect(newUserDefinedGenes).toStrictEqual(geneList.slice(0, globals.maxGenes));
    expect(newDiffExpGenes).toStrictEqual([]);
  });
});
