import React from "react";
import icon from "../../images/icon.png";

const Logo = (props) => {
  const { size } = props;
  return (
    <img
      src={icon}
      height={size}
      width={size}
      alt="CELLxGENE Annotate Logo"
     />
  );
};

export default Logo;
