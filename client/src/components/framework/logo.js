import React from "react";
import * as globals from "../../globals";

const Logo = (props) => {
  const { size } = props;
  return (
    <svg width={size} height={size} viewBox="0 0 48 48" fill="none">
      <rect width="48" height="48" fill="white" />
      <rect x="14" y="14" width="34" height="34" fill={globals.logoColor} />
      <rect x="22" y="22" width="18" height="18" fill="white" />
      <rect x="0" y="14" width="8" height="34" fill={globals.logoColor} />
      <rect x="14" y="0" width="34" height="8" fill={globals.logoIconAccentColor} />
    </svg>
  );
};

export default Logo;
