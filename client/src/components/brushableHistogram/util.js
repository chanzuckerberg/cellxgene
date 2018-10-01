import * as globals from "../../globals";
import _ from "lodash";
import * as constants from "./constants";

export const calcHistogram = (ranges, allValuesForContinuousFieldAsArray) => {
  const histogramData = {};

  histogramData.x = d3
    .scaleLinear()
    .domain([ranges.min, ranges.max])
    .range([0, constants.width]);

  histogramData.y = d3
    .scaleLinear()
    .range([constants.height - constants.marginBottom, 0]);
  // .range([height - margin.bottom, margin.top]);

  histogramData.bins = d3
    .histogram()
    .domain(histogramData.x.domain())
    .thresholds(40)(allValuesForContinuousFieldAsArray);

  histogramData.numValues = allValuesForContinuousFieldAsArray.length;

  console.log(histogramData);
  return histogramData;
};
