import React, { useCallback } from "react";
import { Button, ButtonGroup, Tooltip, Icon } from "@blueprintjs/core";
import { IconNames } from "@blueprintjs/icons";
import * as globals from "../../globals";

const HistogramHeader = React.memo(
  ({
    // @ts-expect-error ts-migrate(2339) FIXME: Property 'fieldId' does not exist on type '{ child... Remove this comment to see the full error message
    fieldId,
    // @ts-expect-error ts-migrate(2339) FIXME: Property 'isColorBy' does not exist on type '{ chi... Remove this comment to see the full error message
    isColorBy,
    // @ts-expect-error ts-migrate(2339) FIXME: Property 'onColorByClick' does not exist on type '... Remove this comment to see the full error message
    onColorByClick,
    // @ts-expect-error ts-migrate(2339) FIXME: Property 'onRemoveClick' does not exist on type '{... Remove this comment to see the full error message
    onRemoveClick,
    // @ts-expect-error ts-migrate(2339) FIXME: Property 'isScatterPlotX' does not exist on type '... Remove this comment to see the full error message
    isScatterPlotX,
    // @ts-expect-error ts-migrate(2339) FIXME: Property 'isScatterPlotY' does not exist on type '... Remove this comment to see the full error message
    isScatterPlotY,
    // @ts-expect-error ts-migrate(2339) FIXME: Property 'onScatterPlotXClick' does not exist on t... Remove this comment to see the full error message
    onScatterPlotXClick,
    // @ts-expect-error ts-migrate(2339) FIXME: Property 'onScatterPlotYClick' does not exist on t... Remove this comment to see the full error message
    onScatterPlotYClick,
    // @ts-expect-error ts-migrate(2339) FIXME: Property 'isObs' does not exist on type '{ childre... Remove this comment to see the full error message
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
            <Icon icon={IconNames.SCATTER_PLOT} style={{ marginRight: 7 }} />
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
