import {
  countCategoryValues2D,
  clearCaches
} from "../../../src/util/stateManager/worldUtil";
import * as Dataframe from "../../../src/util/dataframe";

describe("WorldUtil cache management", () => {
  test("empty", () => {
    const count = countCategoryValues2D(
      "a",
      "b",
      new Dataframe.Dataframe([0, 0], [])
    );
    expect(count).toMatchObject(new Map());
    expect(count.size).toBe(0);
  });

  test("simple couts", () => {
    const df = new Dataframe.Dataframe(
      [3, 2],
      [[0, 0, 1], [false, true, false]],
      null,
      new Dataframe.KeyIndex(["a", "b"])
    );
    const count = countCategoryValues2D("a", "b", df);
    expect(count).toMatchObject(
      new Map([
        [0, new Map([[true, 1], [false, 1]])],
        [1, new Map([[false, 1]])]
      ])
    );
  });

  test("memo cache clear", () => {
    clearCaches();
    const df1 = new Dataframe.Dataframe([0, 0], []);
    const df2 = new Dataframe.Dataframe(
      [3, 2],
      [[0, 0, 1], [false, true, false]],
      null,
      new Dataframe.KeyIndex(["a", "b"])
    );

    const count1 = countCategoryValues2D("a", "b", df1);
    const count2 = countCategoryValues2D("a", "b", df1);
    const count3 = countCategoryValues2D("a", "b", df1.clone());
    const count4 = countCategoryValues2D("a", "b", df2);

    clearCaches();
    const count10 = countCategoryValues2D("a", "b", df1);
    const count11 = countCategoryValues2D("a", "b", df2);

    expect(count1).toEqual(count2);
    expect(count1).toEqual(count3);
    expect(count1).toEqual(count10);
    expect(count1).not.toBe(count3);
    expect(count1).not.toBe(count10);

    expect(count4).toEqual(count11);
    expect(count4).not.toBe(count11);
  });
});
