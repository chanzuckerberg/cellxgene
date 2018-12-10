import {
  countCategoryValues2D,
  clearCaches
} from "../../../src/util/stateManager/worldUtil";

describe("WorldUtil cache management", () => {
  test("empty", () => {
    const count = countCategoryValues2D("a", "b", []);
    expect(count).toMatchObject(new Map());
  });

  test("simple couts", () => {
    const rows = [{ a: 0, b: false }, { a: 0, b: true }, { a: 1, b: false }];
    const count = countCategoryValues2D("a", "b", rows);
    expect(count).toMatchObject(
      new Map([
        [0, new Map([[true, 1], [false, 1]])],
        [1, new Map([[false, 1]])]
      ])
    );
  });

  test("memo cache clear", () => {
    clearCaches();
    const row1 = [];
    const row2 = [{ a: 0, b: false }, { a: 0, b: true }, { a: 1, b: false }];
    const count1 = countCategoryValues2D("a", "b", row1);
    const count2 = countCategoryValues2D("a", "b", row1);
    const count3 = countCategoryValues2D("a", "b", []);
    const count4 = countCategoryValues2D("a", "b", row2);

    clearCaches();
    const count10 = countCategoryValues2D("a", "b", row1);
    const count11 = countCategoryValues2D("a", "b", row2);

    expect(count1).toEqual(count2);
    expect(count1).toEqual(count3);
    expect(count1).toEqual(count10);
    expect(count1).not.toBe(count3);
    expect(count1).not.toBe(count10);

    expect(count4).toEqual(count11);
    expect(count4).not.toBe(count11);
  });
});
