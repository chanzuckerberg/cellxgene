/* core dependencies */
import { Classes } from "@blueprintjs/core";
import React from "react";

function IconAbout() {
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
        d="M7.99999 2.00002C4.68628 2.00002 1.99999 4.68631 1.99999 8.00002C1.99999 11.3137 4.68628 14 7.99999 14C11.3137 14 14 11.3137 14 8.00002C14 4.68631 11.3137 2.00002 7.99999 2.00002ZM0.666656 8.00002C0.666656 3.94993 3.9499 0.666687 7.99999 0.666687C12.0501 0.666687 15.3333 3.94993 15.3333 8.00002C15.3333 12.0501 12.0501 15.3334 7.99999 15.3334C3.9499 15.3334 0.666656 12.0501 0.666656 8.00002ZM7.33332 5.33335C7.33332 4.96516 7.6318 4.66669 7.99999 4.66669H8.00666C8.37485 4.66669 8.67332 4.96516 8.67332 5.33335C8.67332 5.70154 8.37485 6.00002 8.00666 6.00002H7.99999C7.6318 6.00002 7.33332 5.70154 7.33332 5.33335ZM7.99999 7.33335C8.36818 7.33335 8.66666 7.63183 8.66666 8.00002V10.6667C8.66666 11.0349 8.36818 11.3334 7.99999 11.3334C7.6318 11.3334 7.33332 11.0349 7.33332 10.6667V8.00002C7.33332 7.63183 7.6318 7.33335 7.99999 7.33335Z"
        fill="#5C7080"
      />
    </svg>
  );
}
export default IconAbout;