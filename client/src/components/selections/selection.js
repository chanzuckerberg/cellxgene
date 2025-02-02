import React from "react";
import cls from "./selection.css";

export const Selection = ({ emb, name, indexes }) => (
  <div className={cls.root}>
    <span style={{ fontWeight: "bold" }}>{name}({emb})</span>
    <span>{indexes.length} cells</span>
  </div>
);
