import React, { Component } from "react";
import { Tooltip, Position } from "@blueprintjs/core";
import pixelWidth from "string-pixel-width";

import { widthsMap, tooltipHoverOpenDelayQuick } from "../../globals";

export default class Truncate extends Component {
  maybeTruncateString = (str, maxSize, activeFont, fontSize) => {
    let truncatedString = null;
    const renderedSize = pixelWidth(str, {
      font: activeFont,
      size: fontSize,
      map: widthsMap,
    });
    if (renderedSize > maxSize) {
      const percentage = renderedSize / maxSize;
      const newLength = str.length / percentage;
      truncatedString = `${str.slice(0, newLength / 2)}…${str.slice(
        -newLength / 2
      )}`;
    }

    return truncatedString;
  };

  render() {
    const {
      children,
      size,
      fontSize,
      "data-testid": testID,
      "data-testclass": testClass,
      style,
    } = this.props;
    const truncatedString = this.maybeTruncateString(
      children,
      size,
      "Roboto Condensed",
      fontSize
    );
    return (
      <Tooltip
        content={children}
        disabled={truncatedString === null}
        hoverOpenDelay={tooltipHoverOpenDelayQuick}
        position={Position.LEFT}
        usePortal
        modifiers={{
          preventOverflow: { enabled: false },
          hide: { enabled: false },
        }}
      >
        <span data-testid={testID} data-testclass={testClass} style={style}>
          {truncatedString || children}
        </span>
      </Tooltip>
    );
  }
}
