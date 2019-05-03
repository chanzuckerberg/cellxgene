import { range, rangeFill } from "../../src/util/range";

describe("range", () => {
	test("no defaults", () => {
		expect(range(0, 3, 1)).toMatchObject([0, 1, 2]);
	});

	test("range(stop)", () => {
		expect(range(3)).toMatchObject([0, 1, 2]);
		expect(range(0)).toMatchObject([]);
		expect(range(1)).toMatchObject([0]);
	});

	test("range(start,stop)", () => {
		expect(range(0, 0)).toMatchObject([]);
		expect(range(0, 2)).toMatchObject([0, 1]);
		expect(range(4, 8)).toMatchObject([4, 5, 6, 7]);
	});

	test("range(start, stop, step", () => {
		expect(range(4, 0, -1)).toMatchObject([4, 3, 2, 1]);
		expect(range(0, 4, 2)).toMatchObject([0, 2]);
	});
});

describe("rangefill", () => {
	test("rangeFill(arr)", () => {
		expect(rangeFill(new Int32Array(3))).toMatchObject(
			new Int32Array([0, 1, 2])
		);
	});
	test("rangeFill(arr, start)", () => {
		expect(rangeFill(new Int32Array(2), 1)).toMatchObject(
			new Int32Array([1, 2])
		);
	});
	test("rangeFill(arr, start, step)", () => {
		expect(rangeFill(new Int32Array(3), 2, -1)).toMatchObject(
			new Int32Array([2, 1, 0])
		);
	});
});
