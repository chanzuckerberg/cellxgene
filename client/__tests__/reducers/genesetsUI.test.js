import genesetsUIReducer from "../../src/reducers/genesetsUI";

// Format: GeneSetsUI(state,action)

const initialState = {
  createGenesetModeActive: false,
  isEditingGenesetName: false,
  isAddingGenesToGeneset: false,
};

/* initial */
describe("geneset UI states", () => {
  test("initial state, some other action", () => {
    expect(
      genesetsUIReducer(undefined, {
        type: "foo",
      })
    ).toMatchObject(initialState);
  });
  test("geneset: activate add new geneset mode", () => {
    expect(
      genesetsUIReducer(initialState, {
        type: "geneset: activate add new geneset mode",
      })
    ).toMatchObject({
      createGenesetModeActive: true,
      isEditingGenesetName: false,
      isAddingGenesToGeneset: false,
    });
  });
  test("geneset: disable create geneset mode", () => {
    expect(
      genesetsUIReducer(undefined, { isEditingGenesetName: false })
    ).toMatchObject(initialState);
  });

  test("activate add new genes mode", () => {
    expect(
      genesetsUIReducer(undefined, {
        type: "geneset: activate add new genes mode",
        geneset: "a geneset name",
      })
    ).toMatchObject({
      createGenesetModeActive: false,
      isEditingGenesetName: false,
      isAddingGenesToGeneset: "a geneset name",
    });
  });
  test("disable create geneset mode", () => {
    expect(
      genesetsUIReducer(undefined, {
        type: "geneset: disable create geneset mode",
      })
    ).toMatchObject(initialState);
  });
  test("activate rename geneset mode", () => {
    expect(
      genesetsUIReducer(undefined, {
        type: "geneset: activate rename geneset mode",
        data: "a geneset name",
      })
    ).toMatchObject({
      createGenesetModeActive: false,
      isEditingGenesetName: "a geneset name",
      isAddingGenesToGeneset: false,
    });
  });
  test("disable rename geneset mode", () => {
    expect(
      genesetsUIReducer(undefined, {
        type: "geneset: disable rename geneset mode",
      })
    ).toMatchObject(initialState);
  });
});
