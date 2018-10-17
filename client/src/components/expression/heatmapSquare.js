import React from "react";
import getContrast from "font-color-contrast"; // https://www.npmjs.com/package/font-color-contrast

const HeatmapSquare = props => {
  const { backgroundColor, text } = props;
  const contrastColor = getContrast(
    backgroundColor
      .substring(4, backgroundColor.length - 1)
      .replace(/ /g, "")
      .split(",")
  );
  return (
    <p
      style={{
        padding: "12px 6px",
        textAlign: "center",
        color: contrastColor,
        width: 40,
        flexShrink: 0,
        fontSize: 12,
        margin: 0,
        backgroundColor
      }}
    >
      {text}
    </p>
  );
};

export default HeatmapSquare;
