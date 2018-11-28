import WorldOps from "../../../src/util/stateManager/worldOps";

describe("WorldOps cache management", () => {
  test("empty", () => {
    const count = WorldOps.countCategoryValues2D("a", "b", []);
    expect(count).toMatchObject(new Map());
  });

  test("simple couts", () => {
    const rows = [{ a: 0, b: false }, { a: 0, b: true }, { a: 1, b: false }];
    const count = WorldOps.countCategoryValues2D("a", "b", rows);
    expect(count).toMatchObject(
      new Map([
        [0, new Map([[true, 1], [false, 1]])],
        [1, new Map([[false, 1]])]
      ])
    );
  });

  test("memo cache clear", () => {
    WorldOps.clearCaches();
    const row1 = [];
    const row2 = [{ a: 0, b: false }, { a: 0, b: true }, { a: 1, b: false }];
    const count1 = WorldOps.countCategoryValues2D("a", "b", row1);
    const count2 = WorldOps.countCategoryValues2D("a", "b", row1);
    const count3 = WorldOps.countCategoryValues2D("a", "b", []);
    const count4 = WorldOps.countCategoryValues2D("a", "b", row2);

    WorldOps.clearCaches();
    const count10 = WorldOps.countCategoryValues2D("a", "b", row1);
    const count11 = WorldOps.countCategoryValues2D("a", "b", row2);

    expect(count1).toEqual(count2);
    expect(count1).toEqual(count3);
    expect(count1).toEqual(count10);
    expect(count1).not.toBe(count3);
    expect(count1).not.toBe(count10);

    expect(count4).toEqual(count11);
    expect(count4).not.toBe(count11);
  });
});
