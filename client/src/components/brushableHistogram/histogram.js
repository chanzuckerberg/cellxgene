import React, { useEffect, useRef, useState } from "react";
import { interpolateCool } from "d3-scale-chromatic";
import * as d3 from "d3";

import maybeScientific from "./maybeScientific";
import clamp from "./clamp";

const Histogram = ({
  field,
  fieldForId,
  display,
  histogram,
  width,
  height,
  onBrush,
  onBrushEnd,
  margin,
  isColorBy,
  selectionRange,
  mini,
}) => {
  const svgRef = useRef(null);
  const [brush, setBrush] = useState(null);

  useEffect(() => {
    /*
      Create the d3 histogram
      */
    const { marginLeft, marginRight, marginBottom, marginTop } = margin;
    const { x, y, bins, binStart, binEnd, binWidth } = histogram;
    const svg = d3.select(svgRef.current);
    const binPadding = mini ? 0 : -1;
    const defaultBarColor = mini ? "black" : "#bbb";

    /* Remove everything */
    svg.selectAll("*").remove();

    /* Set margins within the SVG */
    const container = svg
      .attr("width", width + marginLeft + marginRight)
      .attr("height", height + marginTop + marginBottom)
      .append("g")
      .attr("class", "histogram-container")
      .attr("transform", `translate(${marginLeft},${marginTop})`);

    const colorScale = d3
      .scaleSequential(interpolateCool)
      .domain([0, bins.length]);

    const histogramScale = d3
      .scaleLinear()
      .domain(x.domain())
      .range([
        colorScale.domain()[1],
        colorScale.domain()[0],
      ]); /* we flip this to make colors dark if high in the color scale */

    if (binWidth > 0) {
      /* BINS */
      container
        .insert("g", "*")
        .selectAll("rect")
        .data(bins)
        .enter()
        .append("rect")
        .attr("x", (d, i) => x(binStart(i)) + 1)
        .attr("y", (d) => y(d))
        .attr("width", (d, i) => x(binEnd(i)) - x(binStart(i)) - binPadding)
        .attr("height", (d) => y(0) - y(d))
        .style(
          "fill",
          isColorBy
            ? (d, i) => colorScale(histogramScale(binStart(i)))
            : defaultBarColor
        );
    }

    // BRUSH
    // Note the brushable area is bounded by the data on three sides, but goes down to cover the x-axis
    const brushX = d3
      .brushX()
      .extent([
        [x.range()[0], y.range()[1]],
        [x.range()[1], marginTop + height + marginBottom],
      ])
      /*
        emit start so that the Undoable history can save an undo point
        upon drag start, and ignore the subsequent intermediate drag events.
        */
      .on("start", onBrush(field, x.invert, "start"))
      .on("brush", onBrush(field, x.invert, "brush"))
      .on("end", onBrushEnd(field, x.invert));

    const brushXselection = container
      .insert("g")
      .attr("class", "brush")
      .attr("data-testid", `${svgRef.current.dataset.testid}-brushable-area`)
      .call(brushX);

    /* X AXIS */
    container
      .insert("g")
      .attr("class", "axis axis--x")
      .attr("transform", `translate(0,${marginTop + height})`)
      .call(
        d3
          .axisBottom(x)
          .ticks(4)
          .tickFormat(d3.format(maybeScientific(x)))
      );

    /* Y AXIS */
    container
      .insert("g")
      .attr("class", "axis axis--y")
      .attr("transform", `translate(${marginLeft + width},0)`)
      .call(
        d3
          .axisRight(y)
          .ticks(3)
          .tickFormat(
            d3.format(
              y.domain().some((n) => Math.abs(n) >= 10000) ? ".0e" : ","
            )
          )
      );

    /* axis style */
    svg.selectAll(".axis text").style("fill", "rgb(80,80,80)");
    svg.selectAll(".axis path").style("stroke", "rgb(230,230,230)");
    svg.selectAll(".axis line").style("stroke", "rgb(230,230,230)");

    setBrush({ brushX, brushXselection });
  }, [histogram, isColorBy]);

  useEffect(() => {
    /*
      paint/update selection brush
      */
    if (!brush) return;
    const { brushX, brushXselection } = brush;
    const selection = d3.brushSelection(brushXselection.node());
    if (!selectionRange && selection) {
      /* no active selection - clear brush */
      brushXselection.call(brushX.move, null);
    } else if (selectionRange) {
      const { x, domain } = histogram;
      const [min, max] = domain;
      const x0 = x(clamp(selectionRange[0], [min, max]));
      const x1 = x(clamp(selectionRange[1], [min, max]));
      if (!selection) {
        /* there is an active selection, but no brush - set the brush */
        brushXselection.call(brushX.move, [x0, x1]);
      } else {
        /* there is an active selection and a brush - make sure they match */
        const moveDeltaThreshold = 1;
        const dX0 = Math.abs(x0 - selection[0]);
        const dX1 = Math.abs(x1 - selection[1]);
        /*
          only update the brush if it is grossly incorrect,
          as defined by the moveDeltaThreshold
          */
        if (dX0 > moveDeltaThreshold || dX1 > moveDeltaThreshold) {
          brushXselection.call(brushX.move, [x0, x1]);
        }
      }
    }
  }, [brush, selectionRange]);

  return (
    <svg
      style={{ display }}
      width={width}
      height={height}
      id={`histogram_${fieldForId}_svg`}
      data-testclass="histogram-plot"
      data-testid={`histogram-${field}-plot`}
      ref={svgRef}
    />
  );
};

export default Histogram;
