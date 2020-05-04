/*
Test cases for nameCreators.js.
*/
import {
  layoutDimensionName,
  obsAnnoDimensionName,
  diffexpDimensionName,
  userDefinedDimensionName,
  makeContinuousDimensionName,
} from "../../src/util/nameCreators";

describe("nameCreators", () => {
  const nameCreators = [
    layoutDimensionName,
    obsAnnoDimensionName,
    diffexpDimensionName,
    userDefinedDimensionName,
  ];

  test("check for namespace isolation", () => {
    const foo = "foo";
    nameCreators.forEach((fn) => expect(fn(foo)).not.toBe(foo));

    const bar = "bar";
    nameCreators.forEach((fn) => expect(fn(foo)).not.toBe(fn(bar)));

    nameCreators.forEach((fn) => {
      const allOtherFn = nameCreators.filter((elmnt) => elmnt !== fn);
      allOtherFn.forEach((otherFn) => expect(fn(foo)).not.toBe(otherFn(foo)));
    });
  });

  test("check for legal keys", () => {
    /* need namespace creators to return strings only */
    nameCreators.forEach((fn) => expect(fn("X")).toMatch(/X/));
  });
});

describe("makeContinuousDimensionName", () => {
  test("namespaces", () => {
    expect(makeContinuousDimensionName({ isObs: true }, "OK")).toMatch(/OK/);
    expect(makeContinuousDimensionName({ isDiffExp: true }, "OK")).toMatch(
      /OK/
    );
    expect(makeContinuousDimensionName({ isUserDefined: true }, "OK")).toMatch(
      /OK/
    );
  });

  test("error handling", () => {
    expect(() => makeContinuousDimensionName({}, "oops")).toThrow();
    expect(() =>
      makeContinuousDimensionName(
        { isObs: false, isDiffExp: false, isUserDefined: false },
        "oops"
      )
    ).toThrow();
    expect(() =>
      makeContinuousDimensionName({ isObs: true }, "OK")
    ).not.toThrow();
    expect(() =>
      makeContinuousDimensionName({ isDiffExp: true }, "OK")
    ).not.toThrow();
    expect(() =>
      makeContinuousDimensionName({ isUserDefined: true }, "OK")
    ).not.toThrow();
  });
});
