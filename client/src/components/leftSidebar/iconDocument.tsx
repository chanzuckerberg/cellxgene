/* core dependencies */
import { Classes } from "@blueprintjs/core";
import React from "react";

// eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types --- FIXME: disabled temporarily on migrate to TS.
function IconDocument() {
  /*
    TODO(cc) Generalize iconography into single component with icon prop.
     */
  return (
    <svg
      className={Classes.ICON}
      fill="none"
      height="16"
      viewBox="0 0 16 16"
      width="16"
      xmlns="http://www.w3.org/2000/svg"
    >
      <path
        fillRule="evenodd"
        clipRule="evenodd"
        d="M2.58579 1.25244C2.96086 0.87737 3.46957 0.666656 4 0.666656H8.66667C8.84348 0.666656 9.01305 0.736894 9.13807 0.861919L13.8047 5.52858C13.9298 5.65361 14 5.82318 14 5.99999V13.3333C14 13.8638 13.7893 14.3725 13.4142 14.7475C13.0391 15.1226 12.5304 15.3333 12 15.3333H4C3.46957 15.3333 2.96086 15.1226 2.58579 14.7475C2.21071 14.3725 2 13.8638 2 13.3333V2.66666C2 2.13622 2.21071 1.62752 2.58579 1.25244ZM4 1.99999C3.82319 1.99999 3.65362 2.07023 3.5286 2.19525C3.40357 2.32028 3.33333 2.48985 3.33333 2.66666V13.3333C3.33333 13.5101 3.40357 13.6797 3.5286 13.8047C3.65362 13.9298 3.82319 14 4 14H12C12.1768 14 12.3464 13.9298 12.4714 13.8047C12.5964 13.6797 12.6667 13.5101 12.6667 13.3333V6.66666H8.66667C8.29848 6.66666 8 6.36818 8 5.99999V1.99999H4ZM9.33333 2.9428L11.7239 5.33332H9.33333V2.9428Z"
        fill="#5C7080"
      />
    </svg>
  );
}
export default IconDocument;
