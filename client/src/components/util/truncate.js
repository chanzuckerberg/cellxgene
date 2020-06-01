import React, { Component, cloneElement } from "react";
import { Tooltip, Position } from "@blueprintjs/core";
import pixelWidth from "string-pixel-width";

import { widthsMap, tooltipHoverOpenDelayQuick } from "../../globals";
import { memoize } from "../../util/dataframe/util";

export default class Truncate extends Component {
  static fontIsLoaded = false;

  // Uses binary-search-esque method to find largest fitting string given a maxSize
  // Params: font, size, bold, italic, map
  // Memoizing requires this function to be above the memo function
  // eslint-disable-next-line react/sort-comp
  static _computeLargestTruncatedString = (
    str,
    maxSize,
    params,
    length = str.length - 1,
    lastLength = 0,
    close = -1
  ) => {
    // Generate the truncated string with ellipses(not if full string)
    const shortenedString =
      length === str.length - 1
        ? str
        : `${str.slice(0, length)}…${str.slice(-length)}`;
    // Measure the size of that string
    const renderedSize = pixelWidth(shortenedString, params);

    // If that full length string is smaller than the maxSize don't return a truncated string
    if (renderedSize < maxSize && length === str.length - 1) {
      return null;
    }
    // Base case: If we've narrowed down to a length, return the closest string
    if (length === lastLength) {
      return `${str.slice(0, close)}…${str.slice(-close)}`;
    }
    // If the current size is longer than max
    if (renderedSize > maxSize) {
      // If we're only one off from the closest string so far, just return the closest string
      if (length - close === 1) {
        return `${str.slice(0, close)}…${str.slice(-close)}`;
      }
      // Otherwise recursively call with half the current length
      return this._computeLargestTruncatedString(
        str,
        maxSize,
        params,
        length < lastLength
          ? Math.floor(length / 2)
          : Math.floor(Math.abs(length - lastLength) / 2) + lastLength,
        length,
        close
      );
    }
    // If the current size is smaller than the max
    if (renderedSize < maxSize) {
      // Save this length as the closest so far
      close = length;
      // Recursively call with a larger string
      return this._computeLargestTruncatedString(
        str,
        maxSize,
        params,
        Math.floor(Math.abs(lastLength - length) / 2) + length,
        length,
        close
      );
    }

    return shortenedString;
  };

  static _memoComputeLargestTruncatedString = memoize(
    this._computeLargestTruncatedString,
    (str, maxSize, params) => {
      return `${str} + ${maxSize} + {${params.font}, ${params.size}, ${params.bold}, ${params.italic}}`;
    }
  );

  // Check's to see if main font (Roboto Condensed) is loaded
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

  // Retrieves the currently loaded font
  static getLoadedFont() {
    if (!this.fontIsLoaded) {
      this.fontIsLoaded = this.checkIfFontLoaded();
    }

    return this.fontIsLoaded ? "Roboto Condensed" : "Helvetica Neue";
  }

  // Returns a truncated string if necessary, null otherwise
  maybeTruncateString = (str, maxSize, fontSize, bold, italic) => {
    const activeFont = Truncate.getLoadedFont();
    return Truncate._memoComputeLargestTruncatedString(str, maxSize, {
      font: activeFont,
      size: fontSize,
      map: widthsMap,
      bold,
      italic,
    });
  };

  render() {
    const { children } = this.props;
    // Truncate only support a single child with a text child
    if (
      React.Children.count(children) === 1 &&
      React.Children.count(children.props?.children) === 1
    ) {
      const originalStr = children.props.children;
      const firstHalf = originalStr.slice(0, originalStr.length / 2);
      const secondHalf = originalStr.slice(originalStr.length / 2);

      const startFixedChars = 5;
      const endFixedChars = 5;
      const fontFaceScaleFactor = 0.16;
      const startWidth = `1em * ${fontFaceScaleFactor} * ${startFixedChars} * 2`;
      const endWidth = `1em * ${fontFaceScaleFactor} * ${endFixedChars} * 2`;

      const commonStyling = {
        display: "inline-block",
        verticalAlign: "bottom",
        whiteSpace: "nowrap",
        overflow: "hidden",
      };

      const truncatedChildren = (
        <>
          <span
            style={{
              ...commonStyling,
              maxWidth: `calc(100% - calc(${endWidth}))`,
              minWidth: `calc(${startWidth})`,
              textOverflow: "ellipsis",
              backgroundColor: "rgba(255, 0, 0, 0.7)",
            }}
          >
            {firstHalf}
          </span>
          <span
            style={{
              ...commonStyling,
              maxWidth: `calc(100% - calc(${startWidth}))`,
              direction: "rtl",
              backgroundColor: "rgba(0, 0, 255, 0.7)",
            }}
          >
            {secondHalf}
          </span>
        </>
      );

      const newChildren = React.Children.map(children, (child) =>
        cloneElement(child, {
          "data-truncated": true,
          children: truncatedChildren,
          style: {
            ...child.props.style,
            whiteSpace: "nowrap",
            overflow: "hidden",
          },
        })
      );
      return (
        <Tooltip
          content={originalStr}
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

    throw Error("Only pass a single child with inner text to Truncate");
  }
}
