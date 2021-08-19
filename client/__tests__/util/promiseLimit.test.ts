import PromiseLimit from "../../src/util/promiseLimit";
import { range } from "../../src/util/range";

const delay = (t: any) => new Promise((resolve) => setTimeout(resolve, t));

describe("PromiseLimit", () => {
  test("simple evaluation, concurrency 1", async () => {
    const plimit = new PromiseLimit(1);
    const result = await Promise.all([
      plimit.add(() => Promise.resolve(1)),
      plimit.add(() => Promise.resolve(2)),
      plimit.add(() => Promise.resolve(3)),
      plimit.add(() => Promise.resolve(4)),
    ]);
    expect(result).toEqual([1, 2, 3, 4]);
  });

  test("simple evaluation, concurrency > 1", async () => {
    const plimit = new PromiseLimit(100);
    const result = await Promise.all([
      plimit.add(() => Promise.resolve(1)),
      plimit.add(() => Promise.resolve(2)),
      plimit.add(() => Promise.resolve(3)),
      plimit.add(() => Promise.resolve(4)),
    ]);
    expect(result).toEqual([1, 2, 3, 4]);
  });

  test("eval in order of insertion", async () => {
    const plimit = new PromiseLimit(100);

    const result = await Promise.all([
      plimit.add(() => Promise.resolve(1)),
      plimit.add(() => Promise.resolve(2)),
      plimit.add(() => Promise.resolve(3)),
      plimit.add(() => Promise.resolve(4)),
    ]);

    expect(result).toEqual([1, 2, 3, 4]);
  });

  test("obeys concurrency limit", async () => {
    const plimit = new PromiseLimit(2);
    let running = 0;
    let maxRunning = 0;

    const callback = async () => {
      running += 1;
      maxRunning = running > maxRunning ? running : maxRunning;
      await delay(100);
      running -= 1;
    };

    // @ts-expect-error ts-migrate(2554) FIXME: Expected 3 arguments, but got 1.
    await Promise.all(range(10).map((i: any) => plimit.add(() => callback(i))));

    expect(maxRunning).toEqual(2);
  });

  test("rejection", async () => {
    const plimit = new PromiseLimit(2);
    const result = await Promise.all([
      plimit.add(() => Promise.resolve("OK")),
      // eslint-disable-next-line prefer-promise-reject-errors -- unit test
      plimit.add(() => Promise.reject("not OK")).catch((e) => e),
      plimit.add(() => Promise.resolve("OK")),
      plimit
        .add(() => {
          throw new Error("not OK");
        })
        .catch((e) => e.message),
    ]);
    expect(result).toEqual(["OK", "not OK", "OK", "not OK"]);
  });

  test("priority queue", async () => {
    const plimit = new PromiseLimit(1);

    let finishOrder = 0;
    const callback = () => async () => {
      await delay(100);
      const result = finishOrder;
      finishOrder += 1;
      return result;
    };

    const result = await Promise.all([
      plimit.add(callback()),
      plimit.priorityAdd(4, callback()),
      plimit.priorityAdd(0, callback()),
      plimit.priorityAdd(1, callback()),
      plimit.priorityAdd(-1, callback()),
    ]);

    expect(result).toEqual([0, 4, 2, 3, 1]);
  });
});
