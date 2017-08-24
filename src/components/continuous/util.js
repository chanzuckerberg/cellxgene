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
    defaultScale: d3.scaleSqrt().range([innerHeight, 0])
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
    key: "Total_reads",
    type: types["Number"],
    scale: d3.scaleSqrt().range([innerHeight, 0])
  },
  {
    key: "Unique_reads",
    type: types["Number"],
    scale: d3.scaleSqrt().range([innerHeight, 0])
  },
  {
    key: "Unique_reads_percent",
    type: types["Number"],
    scale: d3.scaleSqrt().range([innerHeight, 0])
  },
  {
    key: "Genes_detected",
    type: types["Number"],
    scale: d3.scaleSqrt().range([innerHeight, 0])
  },
  {
    key: "ERCC_reads",
    type: types["Number"],
    scale: d3.scaleSqrt().range([innerHeight, 0])
  },
  {
    key: "Non_ERCC_reads",
    type: types["Number"],
    scale: d3.scaleSqrt().range([innerHeight, 0])
  },
  {
    key: "ERCC_to_non_ERCC",
    type: types["Number"],
    scale: d3.scaleSqrt().range([innerHeight, 0])
  },
  {
    key: "Multimapping_reads_percent",
    type: types["Number"],
    scale: d3.scaleSqrt().range([innerHeight, 0])
  },
  {
    key: "Splice_sites_AT.AC",
    type: types["Number"],
    scale: d3.scaleSqrt().range([innerHeight, 0])
  },
  {
    key: "Splice_sites_Annotated",
    type: types["Number"],
    scale: d3.scaleSqrt().range([innerHeight, 0])
  },
  {
    key: "Splice_sites_GC.AG",
    type: types["Number"],
    scale: d3.scaleSqrt().range([innerHeight, 0])
  },
  {
    key: "Splice_sites_GT.AG",
    type: types["Number"],
    scale: d3.scaleSqrt().range([innerHeight, 0])
  },
  {
    key: "Splice_sites_non_canonical",
    type: types["Number"],
    scale: d3.scaleSqrt().range([innerHeight, 0])
  },
  {
    key: "Splice_sites_total",
    type: types["Number"],
    scale: d3.scaleSqrt().range([innerHeight, 0])
  },
  {
    key: "Unmapped_mismatch",
    type: types["Number"],
    scale: d3.scaleSqrt().range([innerHeight, 0])
  },
  {
    key: "Unmapped_other",
    type: types["Number"],
    scale: d3.scaleSqrt().range([innerHeight, 0])
  },
  {
    key: "Unmapped_short",
    type: types["Number"],
    scale: d3.scaleSqrt().range([innerHeight, 0])
  }
];

export const xscale = d3.scalePoint()
    .domain(d3.range(dimensions.length))
    .range([0, width]);

export const yAxis = d3.axisLeft();
