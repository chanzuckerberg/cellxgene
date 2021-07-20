import React from "react";

function Container(props: any) {
  const { children } = props;
  return (
    <div
      className="container"
      style={{
        height: "calc(100vh - (100vh - 100%))",
        width: "calc(100vw - (100vw - 100%))",
        position: "absolute",
        top: 0,
        left: 0,
      }}
    >
      {children}
    </div>
  );
}

export default Container;
