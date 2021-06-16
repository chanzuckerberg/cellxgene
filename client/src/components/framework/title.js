import React from "react";
import * as globals from "../../globals";

const Title = () => {
  return (
    <span
      style={{
        fontSize: 24,
        position: "relative",
        top: -6,
        fontWeight: "bold",
        marginLeft: 5,
        color: globals.logoColor,
        userSelect: "none",
      }}
    >
      cell
      <span
        style={{
          position: "relative",
          top: 1,
          fontWeight: 300,
          fontSize: 24,
        }}
      >
        ×
      </span>
      gene
    </span>
  );
};

export default Title;
