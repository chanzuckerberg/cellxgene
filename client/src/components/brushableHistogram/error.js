import React from "react";
import * as globals from "../../globals";

const ErrorLoading = ({ displayName, zebra }) => (
    <div
      style={{
        backgroundColor: zebra ? globals.lightestGrey : "white",
        fontStyle: "italic",
      }}
    >
      <span>{`Failure loading ${displayName}`}</span>
    </div>
  );

export default ErrorLoading;
