import React from "react";
import * as globals from "../../globals";

const ErrorLoading = ({ error, width, height }) => {
  console.log(error); // log to console as this is an unepected error
  return (
    <div
      style={{
        position: "fixed",
        fontWeight: 500,
        top: height / 2,
        left: globals.leftSidebarWidth + width / 2 - 50,
      }}
    >
      <span>Failure loading dotplot</span>
    </div>
  );
};

export default ErrorLoading;
