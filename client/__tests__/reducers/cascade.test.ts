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
          // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
          currentState: any,
          // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
          action: any,
          // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
          nextSharedState: any,
          // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
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
          // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
          currentState: any,
          // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
          action: any,
          // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
          nextSharedState: any,
          // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
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
