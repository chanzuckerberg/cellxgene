import React, { useState } from "react";
import cls from "./selection.css";
import { Button } from "@blueprintjs/core";
import { Input, Tooltip } from 'antd'

export const Selection = ({ emb, name, indexes, onRemove, onUpdate }) => {
  const [readOnly, setReadOnly] = useState(true)

  const onBlur = (e) => {
    setReadOnly(true)

    const newName = e.target.value;
    if (newName === name) {
      return;
    }

    onUpdate({ name: newName });
  }

  return (
    <div className={cls.root}>
      <Tooltip title={name}>
        <Input size="small" variant="borderless" readOnly={readOnly} onClick={() => setReadOnly(false)} defaultValue={name} onBlur={onBlur} suffix={`(${emb}) ${indexes.length} cells`} />
      </Tooltip>
      <Button small icon="trash" onClick={onRemove} />
    </div>
  );
}
