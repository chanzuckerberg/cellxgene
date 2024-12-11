import genesetsUIReducer from "../../src/reducers/genesetsUI";

// Format: GeneSetsUI(state,action)

const initialState = {
  createGenesetModeActive: false,
  isEditingGenesetName: false,
  sampleset: false,
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
  test("sampleset: activate add new sampleset mode", () => {
    expect(
      genesetsUIReducer(initialState, {
        type: "sampleset: activate add new sampleset mode",
      })
    ).toMatchObject({
      createGenesetModeActive: true,
      isEditingGenesetName: false,
      sampleset: false,
    });
  });
  test("sampleset: disable create sampleset mode", () => {
    expect(
      genesetsUIReducer(undefined, { isEditingGenesetName: false })
    ).toMatchObject(initialState);
  });

  test("activate add new samples mode", () => {
    expect(
      genesetsUIReducer(undefined, {
        type: "sampleset: activate add new samples mode",
        sampleset: "a sampleset name",
      })
    ).toMatchObject({
      createGenesetModeActive: false,
      isEditingGenesetName: false,
      sampleset: "a sampleset name",
    });
  });
  test("disable create sampleset mode", () => {
    expect(
      genesetsUIReducer(undefined, {
        type: "sampleset: disable create sampleset mode",
      })
    ).toMatchObject(initialState);
  });
  test("activate rename sampleset mode", () => {
    expect(
      genesetsUIReducer(undefined, {
        type: "sampleset: activate rename sampleset mode",
        data: "a sampleset name",
      })
    ).toMatchObject({
      createGenesetModeActive: false,
      isEditingGenesetName: "a sampleset name",
      sampleset: false,
    });
  });
  test("disable rename sampleset mode", () => {
    expect(
      genesetsUIReducer(undefined, {
        type: "sampleset: disable rename sampleset mode",
      })
    ).toMatchObject(initialState);
  });
});
