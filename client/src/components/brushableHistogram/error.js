import React from "react";
import * as globals from "../../globals";

const ErrorLoading = ({ displayName, error, zebra }) => {
  console.log(error); // log to console as this is unexpected
  return (
    <div
      style={{
        padding: globals.leftSidebarSectionPadding,
        backgroundColor: zebra ? globals.lightestGrey : "white",
      }}
    >
      <span>{`Failure loading ${displayName}`}</span>
    </div>
  );
};

export default ErrorLoading;
