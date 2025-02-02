import React from "react";
import cls from "./selection.css";
import { Button } from "@blueprintjs/core";

export const Selection = ({ emb, name, indexes, onRemove }) => (
  <div className={cls.root}>
    <span style={{ fontWeight: "bold" }}>{name}({emb})</span>
    <span>{indexes.length} cells</span>
    <Button small icon="trash" onClick={onRemove} />
  </div>
);
