import React from "react";

const StillLoading = ({ width, height }) => {
  /*
      Render a busy/loading indicator
      */
  return (
    <div
      style={{
        position: "fixed",
        fontWeight: 500,
        top: height / 2,
        width,
      }}
    >
      <div
        style={{
          display: "flex",
          justifyContent: "center",
          justifyItems: "center",
          alignItems: "center",
        }}
      >
        <span style={{ fontStyle: "italic" }}>Loading dotplot</span>
      </div>
    </div>
  );
};

export default StillLoading;
