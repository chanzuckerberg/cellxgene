import React from "react";
import { interpolateCool } from "d3-scale-chromatic";
import * as d3 from "d3";

const Dot = (props) => {
  const {
    categoryValue,
    _categoryValueIndex,
    histogramMap,
    _geneSymbol,
    _geneIndex,
    rowColumnSize,
    // dotColorScale,
  } = props;

  const bins = histogramMap.has(categoryValue)
    ? histogramMap.get(categoryValue)
    : new Array(100).fill(0);

  const totalCells = bins.reduce(
    (acc, current) => acc + current
  ); /* SUM REMAINING ELEMENTS */
  /* 
    TODO(colinmegill) #632 this is a heuristic —
    we need to figure out what non expressing means
    and use it, rather than shifting off the first bin
    for prototyping
  */
  bins.shift(); /* MUTATES, REMOVES FIRST ELEMENT */

  const expressing = bins.reduce(
    (acc, current) => acc + current
  ); /* SUM REMAINING ELEMENTS */

  /* TODO(colinmegill) #632 scale between correct dimensions */
  const paddingEquivalentToRowColumnIndexOffset = 8;
  const dotscale = d3
    .scaleLinear()
    .domain([0, 1])
    .range([0, rowColumnSize - paddingEquivalentToRowColumnIndexOffset]);

  const _radius = dotscale(expressing / totalCells);

  return (
    <g
      id={`row_${categoryValue}_${_geneSymbol}`}
      key={`${_categoryValueIndex}_${categoryValue}`}
      transform={`translate(${_geneIndex * rowColumnSize}, ${
        _categoryValueIndex * rowColumnSize
      })`}
    >
      {_geneIndex === 0 && (
        <text
          textAnchor="end"
          style={{ fill: "black", font: "12px Roboto Condensed" }}
        >
          {categoryValue}
        </text>
      )}
      <circle
        r={_radius}
        cx="11"
        cy="-3.5"
        style={{
          fill: interpolateCool(Math.random()),
          fillOpacity: 1,
          stroke: "none",
        }}
      />
    </g>
  );
};

export default Dot;
