import cascadeReducers from "../../src/reducers/cascade";

describe("create", () => {
  test("from Array", () => {
    expect(cascadeReducers([["foo", () => 0]])).toBeInstanceOf(Function);
  });

  test("from Map", () => {
    expect(cascadeReducers(new Map([["foo", () => 0]]))).toBeInstanceOf(
      Function
    );
  });
});

describe("cascade", () => {
  test("expected arguments provided & cascade ordering", () => {
    const topLevelState = {};
    const topLevelAction = { type: "test" };

    const reducer = cascadeReducers([
      [
        "foo",
        (
          currentState: any,
          action: any,
          nextSharedState: any,
          prevSharedState: any
        ) => {
          expect(currentState).toBeUndefined();
          expect(action).toEqual(topLevelAction);
          expect(nextSharedState).toStrictEqual({});
          expect(prevSharedState).toBe(topLevelState);
          return 0;
        },
      ],
      [
        "bar",
        (
          currentState: any,
          action: any,
          nextSharedState: any,
          prevSharedState: any
        ) => {
          expect(currentState).toBeUndefined();
          expect(action).toEqual(topLevelAction);
          expect(nextSharedState).toStrictEqual({ foo: 0 });
          expect(prevSharedState).toBe(topLevelState);
          return 99;
        },
      ],
    ]);

    const nextState = reducer(topLevelState, topLevelAction);
    expect(nextState).toStrictEqual({ foo: 0, bar: 99 });
    expect(topLevelState).toStrictEqual({});
    expect(topLevelAction).toStrictEqual({ type: "test" });
  });
});
