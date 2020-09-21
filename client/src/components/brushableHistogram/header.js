import React, { useCallback } from "react";
import { Button, ButtonGroup, Tooltip } from "@blueprintjs/core";
import * as globals from "../../globals";

const HistogramHeader = React.memo(
  ({
    fieldId,
    isColorBy,
    onColorByClick,
    onRemoveClick,
    isScatterPlotX,
    isScatterPlotY,
    onScatterPlotXClick,
    onScatterPlotYClick,
    isObs,
  }) => {
    /*
        Render the toolbar for the histogram.  Props:
          * fieldId - field identifier, used for various IDs
          * isColorBy - true/false, is this the current color-by
          * onColorByClick - color-by click handler
          * onRemoveClick - optional handler for remove.  Button will not render if not defined.
          * isScatterPlotX - optional, true/false if currently the X scatterplot field
          * isScatterPlotY - optional, true/false if currently the Y scatterplot field
          * onScatterPlotXClick - optional, handler for scatterPlot X button.
          * onScatterPlotYClick - optional, handler for scatterPlot X button.
  
        Scatterplot controls will not render if either handler unspecified.
      */

    const memoizedColorByCallback = useCallback(
      () => onColorByClick(fieldId, isObs),
      [fieldId, isObs]
    );

    return (
      <div
        style={{
          display: "flex",
          justifyContent: "flex-end",
          paddingBottom: "8px",
        }}
      >
        {onScatterPlotXClick && onScatterPlotYClick ? (
          <span>
            <span
              style={{ marginRight: 7 }}
              className="bp3-icon-standard bp3-icon-scatter-plot"
            />
            <ButtonGroup style={{ marginRight: 7 }}>
              <Button
                data-testid={`plot-x-${fieldId}`}
                onClick={onScatterPlotXClick}
                active={isScatterPlotX}
                intent={isScatterPlotX ? "primary" : "none"}
              >
                plot x
              </Button>
              <Button
                data-testid={`plot-y-${fieldId}`}
                onClick={onScatterPlotYClick}
                active={isScatterPlotY}
                intent={isScatterPlotY ? "primary" : "none"}
              >
                plot y
              </Button>
            </ButtonGroup>
          </span>
        ) : null}
        {onRemoveClick ? (
          <Button
            minimal
            onClick={onRemoveClick}
            style={{
              color: globals.blue,
              cursor: "pointer",
              marginLeft: 7,
            }}
          >
            remove
          </Button>
        ) : null}
        <Tooltip
          content="Use as color scale"
          position="bottom"
          hoverOpenDelay={globals.tooltipHoverOpenDelay}
        >
          <Button
            onClick={memoizedColorByCallback}
            active={isColorBy}
            intent={isColorBy ? "primary" : "none"}
            data-testclass="colorby"
            data-testid={`colorby-${fieldId}`}
            icon="tint"
          />
        </Tooltip>
      </div>
    );
  }
);

export default HistogramHeader;
