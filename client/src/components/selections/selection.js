import React from "react";
import cls from "./selection.css";

export const Selection = ({ name, indexes, onMouseOver, onMouseLeave }) => (
  <div
    className={cls.root}
    onFocus={onMouseOver}
    onMouseOver={onMouseOver}
    onMouseLeave={onMouseLeave}
  >
    <span style={{ fontWeight: "bold" }}>{name}</span>
    <span>{indexes.length} cells</span>
  </div>
);
