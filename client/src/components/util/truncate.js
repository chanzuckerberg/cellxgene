import React, { cloneElement } from "react";
import { Tooltip, Position } from "@blueprintjs/core";

import { tooltipHoverOpenDelayQuick } from "../../globals";

export default (props) => {
  const { children } = props;
  // Truncate only support a single child with a text child

  if (
    React.Children.count(children) === 1 &&
    React.Children.count(children.props?.children) === 1
  ) {
    const originalString = children.props.children;

    const firstString = originalString.substr(0, originalString.length / 2);
    const secondString = originalString.substr(originalString.length / 2);

    const splitStyle = {
      ...children.props.style,
      display: "flex",
      overflow: "hidden",
      justifyContent: "flex-start",
    };

    const firstStyle = {
      overflow: "hidden",
      textOverflow: "ellipsis",
      whiteSpace: "nowrap",
      flexShrink: 1,
      minWidth: 5,
    };

    const secondStyleInner = {
      position: "absolute",
      right: 0,
      color: "initial",
    };
    const secondStyle = {
      color: "transparent",
      position: "relative",
      overflow: "hidden",
      whiteSpace: "nowrap",
    };

    const truncatedJSX = (
      <span style={splitStyle}>
        <span style={firstStyle}>{firstString}</span>
        <span style={secondStyle}>
          {secondString}
          <span style={secondStyleInner}>{secondString}</span>
        </span>
      </span>
    );

    // clone children, changing the children(text) to the truncated string
    const newChildren = React.Children.map(children, (child) =>
      cloneElement(child, {
        children: truncatedJSX,
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
  throw Error("Only pass a single child with inner text to Truncate");
};
