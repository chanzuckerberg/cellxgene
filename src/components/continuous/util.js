export const margin = {top: 66, right: 110, bottom: 20, left: 188};
export const width = document.body.clientWidth - margin.left - margin.right;
export const height = 340 - margin.top - margin.bottom;
export const innerHeight = height - 2;

export const devicePixelRatio = window.devicePixelRatio || 1;

export const color = d3.scaleOrdinal()
  .range(["#5DA5B3","#D58323","#DD6CA7","#54AF52","#8C92E8","#E15E5A","#725D82","#776327","#50AB84","#954D56","#AB9C27","#517C3F","#9D5130","#357468","#5E9ACF","#C47DCB","#7D9E33","#DB7F85","#BA89AD","#4C6C86","#B59248","#D8597D","#944F7E","#D67D4B","#8F86C2"]);

export const types = {
  "Number": {
    key: "Number",
    coerce: function(d) { return +d; },
    extent: d3.extent,
    within: function(d, extent, dim) { return extent[0] <= dim.scale(d) && dim.scale(d) <= extent[1]; },
    defaultScale: d3.scaleLinear().range([innerHeight, 0])
  },
  "String": {
    key: "String",
    coerce: String,
    extent: function (data) { return data.sort(); },
    within: function(d, extent, dim) { return extent[0] <= dim.scale(d) && dim.scale(d) <= extent[1]; },
    defaultScale: d3.scalePoint().range([0, innerHeight])
  },
  "Date": {
    key: "Date",
    coerce: function(d) { return new Date(d); },
    extent: d3.extent,
    within: function(d, extent, dim) { return extent[0] <= dim.scale(d) && dim.scale(d) <= extent[1]; },
    defaultScale: d3.scaleTime().range([0, innerHeight])
  }
};

export const dimensions = [
  {
    key: "Total lipid (fat) (g)",
    type: types["Number"],
    scale: d3.scaleSqrt().range([innerHeight, 0])
  },
  {
    key: "Sugars, total (g)",
    type: types["Number"],
    scale: d3.scaleSqrt().range([innerHeight, 0])
  },
  {
    key: "Calcium, Ca (mg)",
    type: types["Number"],
    scale: d3.scaleSqrt().range([innerHeight, 0])
  },
  {
    key: "Sodium, Na (mg)",
    type: types["Number"],
    scale: d3.scaleSqrt().range([innerHeight, 0])
  },
  {
    key: "Phosphorus, P (mg)",
    type: types["Number"],
    scale: d3.scaleSqrt().range([innerHeight, 0])
  },
  {
    key: "Potassium, K (mg)",
    type: types["Number"],
    scale: d3.scaleSqrt().range([innerHeight, 0])
  },
  {
    key: "Thiamin (mg)",
    type: types["Number"],
    scale: d3.scaleSqrt().range([innerHeight, 0])
  },
  {
    key: "Riboflavin (mg)",
    type: types["Number"],
    scale: d3.scaleSqrt().range([innerHeight, 0])
  },
  {
    key: "Niacin (mg)",
    type: types["Number"],
    scale: d3.scaleSqrt().range([innerHeight, 0])
  },
  {
    key: "Iron, Fe (mg)",
    type: types["Number"],
    scale: d3.scaleSqrt().range([innerHeight, 0])
  },
  {
    key: "Magnesium, Mg (mg)",
    type: types["Number"],
    scale: d3.scaleSqrt().range([innerHeight, 0])
  },
  {
    key: "Protein (g)",
    type: types["Number"],
    scale: d3.scaleSqrt().range([innerHeight, 0])
  },
  {
    key: "Zinc, Zn (mg)",
    type: types["Number"],
    scale: d3.scaleSqrt().range([innerHeight, 0])
  },
  {
    key: "Vitamin B-6 (mg)",
    type: types["Number"],
    scale: d3.scaleSqrt().range([innerHeight, 0])
  },
  {
    key: "Vitamin B-12 (mcg)",
    type: types["Number"],
    scale: d3.scaleSqrt().range([innerHeight, 0])
  },
  {
    key: "Folic acid (mcg)",
    type: types["Number"],
    scale: d3.scaleSqrt().range([innerHeight, 0])
  },
  {
    key: "Selenium, Se (mcg)",
    type: types["Number"],
    scale: d3.scaleSqrt().range([innerHeight, 0])
  },
  {
    key: "Vitamin A, IU (IU)",
    type: types["Number"],
    scale: d3.scaleSqrt().range([innerHeight, 0])
  },
  {
    key: "Vitamin K (phylloquinone) (mcg)",
    type: types["Number"],
    scale: d3.scaleSqrt().range([innerHeight, 0])
  },
  {
    key: "Vitamin C, total ascorbic acid (mg)",
    type: types["Number"],
    scale: d3.scaleSqrt().range([innerHeight, 0])
  },
  {
    key: "Vitamin D (IU)",
    type: types["Number"],
    scale: d3.scaleSqrt().range([innerHeight, 0])
  },
  {
    key: "Cholesterol (mg)",
    type: types["Number"],
    scale: d3.scaleSqrt().range([innerHeight, 0])
  },
  {
    key: "Fiber, total dietary (g)",
    type: types["Number"],
    scale: d3.scaleSqrt().range([innerHeight, 0])
  },
  {
    key: "Carbohydrate, by difference (g)",
    type: types["Number"],
    scale: d3.scaleSqrt().range([innerHeight, 0])
  }
];

export const xscale = d3.scalePoint()
    .domain(d3.range(dimensions.length))
    .range([0, width]);

export const yAxis = d3.axisLeft();
