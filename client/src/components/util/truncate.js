import React, { cloneElement } from "react";
import { Tooltip, Position } from "@blueprintjs/core";

import { tooltipHoverOpenDelayQuick } from "../../globals";

const SPLIT_STYLE = {
  display: "flex",
  overflow: "hidden",
  justifyContent: "flex-start",
};

const FIRST_HALF_STYLE = {
  overflow: "hidden",
  textOverflow: "ellipsis",
  whiteSpace: "nowrap",
  flexShrink: 1,
  minWidth: "5px",
};
const SECOND_HALF_STYLE = {
  position: "relative",
  overflow: "hidden",
  whiteSpace: "nowrap",
};
const SECOND_HALF_SPACING_STYLE = {
  color: "transparent",
};

const SECOND_HALF_INNER_STYLE = {
  position: "absolute",
  right: 0,
};

export default (props) => {
  const { children } = props;
  // Truncate only support a single child with a text child

  if (
    React.Children.count(children) !== 1 ||
    React.Children.count(children.props?.children) !== 1
  ) {
    throw Error("Only pass a single child with text to Truncate");
  }
  const originalString = children.props.children.toString();

  let firstString;
  let secondString;

  if (originalString.length === 1) {
    firstString = originalString;
  } else {
    firstString = originalString.substr(0, originalString.length / 2);
    secondString = originalString.substr(originalString.length / 2);
    if (firstString.charAt(firstString.length - 1) === " ") {
      firstString = `${firstString.substr(0, firstString.length - 1)}\u00a0`;
    }
    if (secondString.charAt(0) === " ") {
      secondString = `\u00a0${secondString.substr(1)}`;
    }
  }

  const inheritedColor = children.props.style?.color;

  const splitStyle = { ...children.props.style, ...SPLIT_STYLE, width: "100%" };
  const secondHalfContentStyle = {
    ...SECOND_HALF_INNER_STYLE,
    color: inheritedColor || "inherit",
  };

  const truncatedJSX = (
    <span style={splitStyle}>
      <span style={FIRST_HALF_STYLE}>{firstString}</span>
      <span style={SECOND_HALF_STYLE}>
        <span style={SECOND_HALF_SPACING_STYLE}>{secondString}</span>
        <span style={secondHalfContentStyle}>{secondString}</span>
      </span>
    </span>
  );

  // clone children, changing the children(text) to the truncated string
  const newChildren = React.Children.map(children, (child) =>
    cloneElement(child, {
      children: truncatedJSX,
      "aria-label": originalString,
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
      targetProps={{ style: children.props.style }}
    >
      {newChildren}
    </Tooltip>
  );
};
