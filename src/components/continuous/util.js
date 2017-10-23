import _ from "lodash";

const paddingRight = 120;

export const margin = {top: 66, right: 110, bottom: 20, left: 60};
export const width = document.body.clientWidth - margin.left - margin.right - paddingRight;
export const height = 340 - margin.top - margin.bottom;
export const innerHeight = height - 2;

export const devicePixelRatio = window.devicePixelRatio || 1;

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

export const createDimensions = (data) => {
  const newArr = []
  _.each(data, (value, key) => {
    if (value.range) {
      newArr.push({
        key: key, /* room for confusion: lodash calls this key, it's also the name of the property parallel coords code is looking for */
        type: types["Number"],
        scale: d3.scaleSqrt().range([innerHeight, 0])
      })
    }
  })
  return newArr;
}

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
