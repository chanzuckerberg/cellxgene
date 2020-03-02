/*
test controls helpers
*/
import { subsetAndResetGeneLists } from "../../../src/util/stateManager/controlsHelpers";
import * as globals from "../../../src/globals";

describe("controls helpers", () => {

  test("subsetAndResetGeneLists", () => {
    const geneList = [
      "Ak7",
      "Agt",
      "Alb",
      "Aqp6",
      "Aqp8",
      "Bves",
      "Bhmt",
      "Bcan",
      "Bfar",
      "Bpgm",
      "Cfd",
      "Clu",
      "Cel",
      "Ckm",
      "Cnp",
      "Ddo",
      "Dld",
      "Ddn",
      "Des",
      "Dmd",
      "Egf",
      "Epc1",
      "Expi",
      "Eno2",
      "Eno3",
      "Fga",
      "Fos",
      "Fgb",
      "Fgg",
      "Fmn2",
    ];
    const state = {
      userDefinedGenes: geneList.slice(0, 15),
      diffexpGenes: geneList.slice(15),
    };
    const [newUserDefinedGenes, newDiffExpGenes] = subsetAndResetGeneLists(state);
    expect(geneList.length).toBeGreaterThan(globals.maxUserDefinedGenes);
    expect(newUserDefinedGenes).toHaveLength(globals.maxUserDefinedGenes);
    expect(newUserDefinedGenes).toStrictEqual(geneList.slice(0, globals.maxUserDefinedGenes));
    expect(newDiffExpGenes).toStrictEqual([]);
  });
});
