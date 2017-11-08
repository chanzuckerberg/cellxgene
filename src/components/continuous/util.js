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

export const yAxis = d3.axisLeft();

export const brushstart = () => {
  d3.event.sourceEvent.stopPropagation();
}

export const d3_functor = (v) => {
  return typeof v === "function" ? v : () => { return v; };
};

export const project = (d, dimensions, xscale) => {
  return dimensions.map((p,i) => {
    // check if data element has property and contains a value
    if (
      !(p.key in d) ||
      d[p.key] === null
    ) return null;

    return [xscale(i),p.scale(d[p.key])];
  });
};

/*****************************************
******************************************
    process data and extract features
******************************************
******************************************/

export const processData = (metadata, dimensions) => {
  // shuffle the metadata! - or not, this is a visual effect
  metadata = d3.shuffle(metadata);

  metadata.forEach((d) => {
    dimensions.forEach((p) => {
      d[p.key] = !d[p.key] ? null : p.type.coerce(d[p.key]);
    });

    // truncate long text strings to fit in metadata table
    for (var key in d) {
      if (d[key] && d[key].length > 35) d[key] = d[key].slice(0,36);
    }
  });

  // type/dimension default setting happens here
  dimensions.forEach((dim) => {
    if (!("domain" in dim)) {
      // detect domain using dimension type's extent function
      dim.domain = d3_functor(dim.type.extent)(
        metadata.map((d) => { return d[dim.key]; })
      );
    }
    if (!("scale" in dim)) {
      // use type's default scale for dimension
      dim.scale = dim.type.defaultScale.copy();
    }
    dim.scale.domain(dim.domain);
  });

  return {
    processedMetadata: metadata,
    processedDimensions: dimensions
  }
}
