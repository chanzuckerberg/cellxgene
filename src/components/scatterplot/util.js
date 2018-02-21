import _ from "lodash";

const paddingRight = 120;
const continuousChartWidth = 1200;

export const margin = {top: 66, right: 110, bottom: 20, left: 60};
export const width = continuousChartWidth - margin.left - margin.right - paddingRight;
export const height = 680 - margin.top - margin.bottom;
export const innerHeight = height - 2;

export const devicePixelRatio = window.devicePixelRatio || 1;
