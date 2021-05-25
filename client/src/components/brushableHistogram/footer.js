import React from "react";

const HistogramFooter = React.memo(
  ({
    displayName,
    hideRanges,
    rangeMin,
    rangeMax,
    rangeColorMin,
    rangeColorMax,
    isObs,
    isGeneSetSummary,
  }) => {
    /*
    Footer of each histogram.  Will render range and title.

    Required props:
      * displayName - the displayName, aka "n_genes", "FOXP2", etc.
      * hideRanges - true/false, enables/disable rendering of ranges
      * range - length two array, [min, max], containing the range values to display
      * rangeColor - length two array, [mincolor, maxcolor], each a CSS color
    */
    return (
      <div>
        <div
          style={{
            display: "flex",
            justifyContent: hideRanges ? "center" : "space-between",
          }}
        >
          <span
            style={{
              color: rangeColorMin,
              display: hideRanges ? "none" : "block",
            }}
          >
            min {rangeMin.toPrecision(4)}
          </span>
          <span
            data-testclass="brushable-histogram-field-name"
            style={{ fontStyle: "italic" }}
          >
            {isObs && displayName}
            {isGeneSetSummary && "gene set mean expression"}
          </span>
          <div style={{ display: hideRanges ? "block" : "none" }}>
            : {rangeMin}
          </div>
          <span
            style={{
              color: rangeColorMax,
              display: hideRanges ? "none" : "block",
            }}
          >
            max {rangeMax.toPrecision(4)}
          </span>
        </div>
      </div>
    );
  }
);

export default HistogramFooter;
