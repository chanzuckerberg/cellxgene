import React from "react";
import {
  Button,
  ButtonGroup,
  Icon,
  Intent,
  NumericInput,
  Popover,
  Position,
  Tooltip,
} from "@blueprintjs/core";
import { IconNames } from "@blueprintjs/icons";

import { tooltipHoverOpenDelay } from "../../globals";
import styles from "./menubar.css";

const Clip = React.memo((props) => {
  const {
    pendingClipPercentiles,
    clipPercentileMin,
    clipPercentileMax,
    handleClipOpening,
    handleClipClosing,
    handleClipCommit,
    isClipDisabled,
    handleClipOnKeyPress,
    handleClipPercentileMaxValueChange,
    handleClipPercentileMinValueChange,
  } = props;

  const clipMin =
    pendingClipPercentiles?.clipPercentileMin ?? clipPercentileMin;
  const clipMax =
    pendingClipPercentiles?.clipPercentileMax ?? clipPercentileMax;
  const intent =
    clipPercentileMin > 0 || clipPercentileMax < 100
      ? Intent.INTENT_WARNING
      : "";

  return (
    <ButtonGroup className={`${styles.menubarButton}`}>
      <Popover
        target={
          <Tooltip
            content="Clip all continuous values to a percentile range"
            position="bottom"
            hoverOpenDelay={tooltipHoverOpenDelay}
          >
            <Button
              type="button"
              data-testid="visualization-settings"
              intent={intent}
              icon={IconNames.TIMELINE_BAR_CHART}
              style={{
                cursor: "pointer",
              }}
            />
          </Tooltip>
        }
        position={Position.BOTTOM_RIGHT}
        onOpening={handleClipOpening}
        onClosing={handleClipClosing}
        content={
          <div
            style={{
              display: "flex",
              justifyContent: "flex-start",
              alignItems: "flex-start",
              flexDirection: "column",
              padding: 10,
            }}
          >
            <div>Clip all continuous values to percentile range</div>
            <div
              style={{
                display: "flex",
                justifyContent: "space-between",
                alignItems: "center",
                paddingTop: 5,
                paddingBottom: 5,
              }}
            >
              <NumericInput
                style={{ width: 50 }}
                data-testid="clip-min-input"
                onValueChange={handleClipPercentileMinValueChange}
                onKeyPress={handleClipOnKeyPress}
                value={clipMin}
                min={0}
                max={100}
                fill={false}
                minorStepSize={null}
                rightElement={
                  <div style={{ padding: "4px 2px" }}>
                    <Icon icon="percentage" intent="primary" iconSize={14} />
                  </div>
                }
              />
              <span style={{ marginRight: 5, marginLeft: 5 }}> - </span>
              <NumericInput
                style={{ width: 50 }}
                data-testid="clip-max-input"
                onValueChange={handleClipPercentileMaxValueChange}
                onKeyPress={handleClipOnKeyPress}
                value={clipMax}
                min={0}
                max={100}
                fill={false}
                minorStepSize={null}
                rightElement={
                  <div style={{ padding: "4px 2px" }}>
                    <Icon icon="percentage" intent="primary" iconSize={14} />
                  </div>
                }
              />
              <Button
                type="button"
                data-testid="clip-commit"
                intent="primary"
                disabled={isClipDisabled()}
                style={{
                  cursor: "pointer",
                  marginRight: 5,
                  marginLeft: 5,
                }}
                onClick={handleClipCommit}
              >
                Clip
              </Button>
            </div>
          </div>
        }
      />
    </ButtonGroup>
  );
});

export default Clip;
