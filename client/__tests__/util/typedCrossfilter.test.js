// jshint esversion: 6
const _ = require("lodash");
const crossfilter = require("../../src/util/typedCrossfilter");

const someData = [
  {
    date: "2011-11-14T16:17:54Z",
    quantity: 2,
    total: 190,
    tip: 100,
    type: "tab",
    productIDs: ["001"]
  },
  {
    date: "2011-11-14T16:20:19Z",
    quantity: 2,
    total: 190,
    tip: 100,
    type: "tab",
    productIDs: ["001", "005"]
  },
  {
    date: "2011-11-14T16:28:54Z",
    quantity: 1,
    total: 300,
    tip: 200,
    type: "visa",
    productIDs: ["004", "005"]
  },
  {
    date: "2011-11-14T16:30:43Z",
    quantity: 2,
    total: 90,
    tip: 0,
    type: "tab",
    productIDs: ["001", "002"]
  },
  {
    date: "2011-11-14T16:48:46Z",
    quantity: 2,
    total: 90,
    tip: 0,
    type: "tab",
    productIDs: ["005"]
  },
  {
    date: "2011-11-14T16:53:41Z",
    quantity: 2,
    total: 90,
    tip: 0,
    type: "tab",
    productIDs: ["001", "004", "005"]
  },
  {
    date: "2011-11-14T16:54:06Z",
    quantity: 1,
    total: 100,
    tip: 0,
    type: "cash",
    productIDs: ["001", "002", "003", "004", "005"]
  },
  {
    date: "2011-11-14T16:58:03Z",
    quantity: 2,
    total: 90,
    tip: 0,
    type: "tab",
    productIDs: ["001"]
  },
  {
    date: "2011-11-14T17:07:21Z",
    quantity: 2,
    total: 90,
    tip: 0,
    type: "tab",
    productIDs: ["004", "005"]
  },
  {
    date: "2011-11-14T17:22:59Z",
    quantity: 2,
    total: 90,
    tip: 0,
    type: "tab",
    productIDs: ["001", "002", "004", "005"]
  },
  {
    date: "2011-11-14T17:25:45Z",
    quantity: 2,
    total: 200,
    tip: 0,
    type: "cash",
    productIDs: ["002"]
  },
  {
    date: "2011-11-14T17:29:52Z",
    quantity: 1,
    total: 200,
    tip: 100,
    type: "visa",
    productIDs: ["004"]
  }
];

var payments = null;
beforeEach(() => {
  payments = crossfilter(someData);
});

describe("typedCrossfilter", () => {
  test("alloc and free", () => {
    expect(payments).toBeDefined();
    expect(payments.size()).toEqual(someData.length);
    expect(payments.all()).toEqual(someData);

    const quantity = payments.dimension(r => r.quantity, Int32Array);
    expect(quantity).toBeDefined();
    expect(quantity.id()).toBeDefined();

    quantity.dispose();
    expect(payments.size()).toEqual(someData.length);
    expect(payments.all()).toEqual(someData);
  });

  test("filterAll and filterNone", () => {
    expect(payments).toBeDefined();
    const quantity = payments.dimension(r => r.quantity, Int32Array);
    const tip = payments.dimension(r => r.tip, Float32Array);
    const total = payments.dimension(r => r.total, Float32Array);
    const type = payments.dimension(r => r.type, "enum");

    expect(quantity).toBeDefined();
    expect(tip).toBeDefined();
    expect(total).toBeDefined();
    expect(type).toBeDefined();

    // initially, all should be filtered
    expect(payments.allFiltered().length).toEqual(payments.size());
    expect(payments.allFiltered()).toEqual(payments.all());
    expect(payments.countFiltered()).toEqual(someData.length);

    // filterAll
    tip.filterAll(); // should change nothing
    expect(payments.allFiltered()).toEqual(payments.all());
    expect(payments.countFiltered()).toEqual(someData.length);

    // ditto
    total.filterAll();
    expect(payments.allFiltered()).toEqual(payments.all());
    expect(payments.countFiltered()).toEqual(someData.length);

    // filterNone
    type.filterNone();
    expect(payments.allFiltered()).toEqual([]);
    expect(payments.countFiltered()).toEqual(0);

    quantity.filterNone();
    expect(payments.allFiltered()).toEqual([]);
    expect(payments.countFiltered()).toEqual(0);

    // invert the first none; should have no effect because type is
    // still not filtered
    quantity.filterAll();
    expect(payments.allFiltered()).toEqual([]);
    expect(payments.countFiltered()).toEqual(0);

    // filter all of type; should select all
    type.filterAll();
    expect(payments.allFiltered()).toEqual(payments.all());
    expect(payments.countFiltered()).toEqual(payments.size());
  });

  test("filterExact", () => {
    expect(payments).toBeDefined();
    const quantity = payments.dimension(r => r.quantity, Int32Array);
    const tip = payments.dimension(r => r.tip, Float32Array);
    const total = payments.dimension(r => r.total, Float32Array);
    const type = payments.dimension(r => r.type, "enum");

    quantity.filterExact(1);
    expect(payments.countFiltered()).toEqual(
      _.countBy(someData, "quantity")[1]
    );
    expect(payments.allFiltered()).toEqual(_.filter(someData, { quantity: 1 }));

    tip.filterExact(0);
    expect(payments.allFiltered()).toEqual(
      _.filter(someData, { tip: 0, quantity: 1 })
    );

    type.filterExact("cash");
    expect(payments.allFiltered()).toEqual(
      _.filter(someData, { tip: 0, quantity: 1, type: "cash" })
    );
  });

  test("filterRange", () => {
    expect(payments).toBeDefined();
    const quantity = payments.dimension(r => r.quantity, Int32Array);
    const tip = payments.dimension(r => r.tip, Float32Array);
    const total = payments.dimension(r => r.total, Float32Array);
    const type = payments.dimension(r => r.type, "enum");

    tip.filterRange([0, 91]);
    expect(payments.allFiltered()).toEqual(
      _(someData)
        .filter(r => r.tip >= 0 && r.tip < 91)
        .value()
    );

    tip.filterRange([0, 90]);
    expect(payments.allFiltered()).toEqual(
      _(someData)
        .filter(r => r.tip >= 0 && r.tip < 90)
        .value()
    );

    tip.filterRange([1, 90]);
    expect(payments.allFiltered()).toEqual(
      _(someData)
        .filter(r => r.tip >= 1 && r.tip < 91)
        .value()
    );
  });

  test("filterEnum", () => {
    expect(payments).toBeDefined();
    const quantity = payments.dimension(r => r.quantity, Int32Array);
    const tip = payments.dimension(r => r.tip, Float32Array);
    const total = payments.dimension(r => r.total, Float32Array);
    const type = payments.dimension(r => r.type, "enum");

    type.filterEnum(["tab", "cash"]);
    expect(payments.allFiltered()).toEqual(
      _(someData)
        .filter(r => r.type === "cash" || r.type === "tab")
        .value()
    );

    tip.filterEnum([0, 100]);
    expect(payments.allFiltered()).toEqual(
      _(someData)
        .filter(r => r.type === "cash" || r.type === "tab")
        .filter(r => r.tip === 0 || r.tip === 100)
        .value()
    );
  });

  test("more than 32 dimensions", () => {
    expect(payments).toBeDefined();
    const quantity = payments.dimension(r => r.quantity, Int32Array);
    const tip = payments.dimension(r => r.tip, Float32Array);
    const total = payments.dimension(r => r.total, Float32Array);
    const type = payments.dimension(r => r.type, "enum");

    // Create a bunch of fake dimensions to ensure we can handle > 32
    let dimMap = {};
    for (let i = 0; i < 65; i++) {
      dimMap[i] = payments.dimension(r => Math.random(), Float32Array);
      expect(dimMap[i]).toBeDefined();
      expect(dimMap[i].id()).toBeDefined();
    }

    // everything should start as selected/filtered
    expect(payments.countFiltered()).toEqual(someData.length);

    dimMap[0].filterAll();
    dimMap[64].filterAll();
    expect(payments.countFiltered()).toEqual(someData.length);

    dimMap[33].filterNone();
    expect(payments.allFiltered()).toEqual([]);

    dimMap[33].filterAll();
    expect(payments.allFiltered()).toEqual(someData);
  });
});
