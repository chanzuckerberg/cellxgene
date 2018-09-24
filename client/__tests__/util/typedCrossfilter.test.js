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

function groupReduce(data, valueMap, valueReduce, valueInit) {
  return _
    .reduce(
      data,
      (acc, value) => {
        const k = valueMap(value);
        let r = _.find(acc, o => o.key === k);
        if (!r) {
          r = { key: k, value: valueInit() };
          acc.push(r);
        }
        r.value = valueReduce(r.value, value);
        return acc;
      },
      []
    )
    .sort((a, b) => (a.key < b.key ? -1 : a.key > b.key ? 1 : 0));
}

function groupCount(data, map) {
  return groupReduce(data, map, (p, v) => p + 1, () => 0);
}

function groupSum(data, map) {
  return groupReduce(data, map, (p, v) => (p += map(v)), () => 0);
}

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

  test("group, default mapping, default reducer, no filter", () => {
    expect(payments).toBeDefined();

    var quantity = payments.dimension(r => r.quantity, Int32Array);
    var tip = payments.dimension(r => r.tip, Int32Array);
    var type = payments.dimension(r => r.type, "enum");
    var total = payments.dimension(r => r.total, Int32Array);

    _.each(
      {
        tip: tip.group(r => r),
        type: type.group(),
        total: total.group(),
        quantity: quantity.group()
      },
      (grp, k) => {
        const whatWeExpect = groupCount(someData, v => v[k]);
        expect(grp.all()).toEqual(whatWeExpect);
        expect(grp.size()).toEqual(whatWeExpect.length);
        expect(grp.dispose()).toEqual(grp);
      }
    );
  });

  test("group, custom map, default reducer, no filters", () => {
    expect(payments).toBeDefined();

    // custom mapping in groups only works for scalar types.  Enums do not
    // currently implement it.

    const tip = payments.dimension(r => r.tip, Int32Array);
    const totalX10 = payments.dimension(r => r.total * 10, Int32Array);
    const type = payments.dimension(r => r.type, "enum");

    const paymentsByTip_A = tip.group();
    const paymentsByTip_B = tip.group(r => 10 * r);
    const paymentsByType = type.group(); // identity only
    const paymentsByTotalX10_A = totalX10.group();
    const paymentsByTotalX10_B = totalX10.group(r => r / 10);

    expect(paymentsByTip_A.all()).toEqual(groupCount(someData, v => v.tip));
    expect(paymentsByTip_B.all()).toEqual(
      groupCount(someData, v => 10 * v.tip)
    );
    expect(paymentsByType.all()).toEqual(groupCount(someData, v => v.type));
    expect(paymentsByTotalX10_A.all()).toEqual(
      groupCount(someData, v => 10 * v.total)
    );
    expect(paymentsByTotalX10_B.all()).toEqual(
      groupCount(someData, v => (10 * v.total) / 10)
    );

    for (let i of [
      paymentsByTip_A,
      paymentsByTip_B,
      paymentsByType,
      paymentsByTotalX10_A,
      paymentsByTotalX10_B,
      tip,
      totalX10,
      type
    ]) {
      expect(i.dispose()).toEqual(i);
    }
  });

  test("group, default map, custom reducer, no filters", () => {
    expect(payments).toBeDefined();

    const total = payments.dimension(r => r.total, Float32Array);
    const type = payments.dimension(r => r.type, "enum");

    const paymentsByTotal = total.group();
    const paymentsByType = type.group();

    // reduceCount
    expect(paymentsByTotal.reduceCount()).toEqual(paymentsByTotal);
    expect(paymentsByTotal.all()).toEqual(groupCount(someData, v => v.total));

    // reduceSum
    expect(paymentsByTotal.reduceSum(v => v.total)).toEqual(paymentsByTotal);
    expect(paymentsByTotal.all()).toEqual(groupSum(someData, v => v.total));

    // use custom reducers (my reducers) - count by three, init 1
    expect(
      paymentsByTotal.reduce((p, v) => (p += 3), (p, v) => (p -= 3), () => 1)
    ).toEqual(paymentsByTotal);
    expect(paymentsByTotal.all()).toEqual(
      groupReduce(someData, v => v.total, (p, v) => p + 3, () => 1)
    );

    for (let i of [paymentsByTotal, paymentsByType, type]) {
      expect(i.dispose()).toEqual(i);
    }
  });

  test("group, default map, default reducer, filters", () => {
    // From the docs:
    // Note: a grouping intersects the crossfilter's current filters, except for the
    // associated dimension's filter. Thus, group methods consider only records that
    // satisfy every filter except this dimension's filter. So, if the crossfilter of
    // payments is filtered by type and total, then group by total only observes the
    // filter by type.

    expect(payments).toBeDefined();

    const tip = payments.dimension(r => r.tip, Int32Array);
    const total = payments.dimension(r => r.total, Int32Array);
    const type = payments.dimension(r => r.type, "enum");

    const paymentsByTip = tip.group();
    const paymentsByTotal = total.group();
    const paymentsByType = type.group();

    // 1. confirm that changing the filter on a dimension does NOT change that
    // dimensions groups.
    {
      tip.filterAll(), total.filterAll(), type.filterAll();
      let before = _.cloneDeep(paymentsByTip.all());
      tip.filterExact(0);
      expect(paymentsByTip.all()).toEqual(before);
    }

    // 2. confirm that changing a filter on a different dimension DOES change
    // all other groups.
    {
      tip.filterAll(), total.filterAll(), type.filterAll();
      const before = _.cloneDeep([paymentsByTotal.all(), paymentsByType.all()]);
      tip.filterExact(0);
      const after = [paymentsByTotal.all(), paymentsByType.all()];
      expect(after).not.toEqual(before);
      expect(after).toEqual([
        groupReduce(
          someData,
          v => v.total,
          (p, v) => (v.tip !== 0 ? p : p + 1),
          () => 0
        ),
        groupReduce(
          someData,
          v => v.type,
          (p, v) => (v.tip !== 0 ? p : p + 1),
          () => 0
        )
      ]);
    }

    for (let i of [
      paymentsByTip,
      paymentsByTotal,
      paymentsByType,
      tip,
      total,
      type
    ]) {
      expect(i.dispose()).toEqual(i);
    }
  });
});
