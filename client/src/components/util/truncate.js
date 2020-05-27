import React, { Component } from "react";
import { Tooltip, Position } from "@blueprintjs/core";
import pixelWidth from "string-pixel-width";

import { widthsMap, tooltipHoverOpenDelayQuick } from "../../globals";

export default class Truncate extends Component {
  static fontIsLoaded = false;

  static checkIfFontLoaded() {
    if (document.fonts.check("1em Roboto Condensed")) {
      return true;
    }
    let ret;
    document.fonts.ready.then(() => {
      if (!document.fonts.check("1em Roboto Condensed")) {
        console.error("Roboto Condensed was not loaded");
        ret = false;
      } else {
        ret = true;
      }
    });
    return ret;
  }

  static getLoadedFont() {
    if (!this.fontIsLoaded) {
      this.fontIsLoaded = this.checkIfFontLoaded();
    }

    return this.fontIsLoaded ? "Roboto Condensed" : "Helvetica Neue";
  }

  maybeTruncateString = (str, maxSize, fontSize) => {
    const activeFont = Truncate.getLoadedFont();

    const renderedSize = pixelWidth(str, {
      font: activeFont,
      size: fontSize,
      map: widthsMap,
    });
    let truncatedString = null;
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
    const truncatedString = this.maybeTruncateString(children, size, fontSize);
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
