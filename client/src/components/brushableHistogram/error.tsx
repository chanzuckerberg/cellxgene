import React from "react";
import * as globals from "../../globals";

// eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any -- - FIXME: disabled temporarily on migrate to TS.
const ErrorLoading = ({ displayName, zebra }: any) => (
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
