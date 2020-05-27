import React, { Component, cloneElement } from "react";
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
    const { children, size, fontSize } = this.props;
    // Truncate only support a single child with a text child

    if (
      React.Children.count(children) === 1 &&
      React.Children.count(children.props?.children) === 1
    ) {
      const originalString = children.props.children;
      const truncatedString = this.maybeTruncateString(
        originalString,
        size,
        fontSize
      );
      // Only make tooltip if string has to be truncated
      if (truncatedString) {
        // clone children, changing the children(text) to the truncated string
        const newChildren = React.Children.map(children, (child) =>
          cloneElement(child, {
            children: truncatedString,
            "data-truncated": true,
          })
        );
        return (
          <Tooltip
            content={originalString}
            hoverOpenDelay={tooltipHoverOpenDelayQuick}
            position={Position.LEFT}
            usePortal
            modifiers={{
              preventOverflow: { enabled: false },
              hide: { enabled: false },
            }}
          >
            {newChildren}
          </Tooltip>
        );
      }
      return children;
    }
    throw Error("Only pass a single child with inner text to Truncate");
  }
}
