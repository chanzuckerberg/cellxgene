import undoable from "../../src/reducers/undoable";

describe("create", () => {
  test("no keys", () => {
    expect(() => undoable(() => {})).toThrow();
    expect(() => undoable(() => {}, null)).toThrow();
    expect(() => undoable(() => {}, [])).toThrow();
    expect(() => undoable(() => {}, [], {})).toThrow();
  });

  test("simple", () => {
    expect(undoable(() => {}, ["foo"])).toBeInstanceOf(Function);
    expect(undoable(() => {}, ["foo"], {})).toBeInstanceOf(Function);
  });

  test("handles undefined initial state", () => {
    expect(
      undoable(() => {}, ["a"])(undefined, { type: "test" })
    ).toMatchObject({});
  });
});

describe("undo", () => {
  test("expected state modifications", () => {
    const initialState = { a: 0, b: 1000 };
    const reducer = (state) => ({ a: state.a + 1, b: state.b + 1 });
    const undoableReducer = undoable(reducer, ["a"]);

    const s1 = undoableReducer(initialState, { type: "test" });
    expect(s1).toMatchObject({ a: 1, b: 1001 });

    // test that only specified keys are undone
    const s2 = undoableReducer(s1, { type: "@@undoable/undo" });
    expect(s2).toMatchObject({ a: 0, b: 1001 });

    // test backstop when no more history
    const s3 = undoableReducer(s2, { type: "@@undoable/undo" });
    expect(s3).toMatchObject({ a: 0, b: 1001 });
  });
});

describe("redo", () => {
  const initialState = { a: 0, b: 1000 };
  const reducer = (state) => ({ a: state.a + 1, b: state.b + 1 });
  let UR;

  beforeEach(() => {
    UR = undoable(reducer, ["a"]);
  });

  test("expected state modifications", () => {
    const s1 = UR(initialState, { type: "test" });
    expect(s1).toMatchObject({ a: 1, b: 1001 });

    // verify undo->redo reverts state.
    const s2 = UR(UR(s1, { type: "@@undoable/undo" }), {
      type: "@@undoable/redo",
    });
    expect(s2).toMatchObject({ a: 1, b: 1001 });

    // verify backstop when no redo future
    const s3 = UR(s2, { type: "@@undoable/redo" });
    expect(s3).toMatchObject({ a: 1, b: 1001 });
  });

  test("history cleared", () => {
    // verify future cleared upon a normal state transition
    const s1 = UR(initialState, { type: "test" });
    expect(s1).toMatchObject({ a: 1, b: 1001 });
    const s2 = UR(s1, { type: "@@undoable/undo" });
    expect(s2).toMatchObject({ a: 0, b: 1001 });
    const s3 = UR(s2, { type: "test" });
    expect(s3).toMatchObject({ a: 1, b: 1002 });
    const s4 = UR(s3, { type: "@@undoable/redo" });
    expect(s4).toMatchObject({ a: 1, b: 1002 });
  });
});

/*
TODO:
- historyLimit is enforced
- action filters
*/
