import React, { useState } from 'react'
import { Button } from '@blueprintjs/core'
import { Input, Tooltip } from 'antd'

export const SelectionGroupTitle = ({ name, onRemove, onUpdate }) => {
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
    <div style={{ display: 'flex', alignItems: 'center', justifyContent: 'space-between', gap: '1ch' }}>
      <Tooltip title={name}>
        <Input size="small" variant="borderless" readOnly={readOnly} onClick={() => setReadOnly(false)} defaultValue={name} onBlur={onBlur} />
      </Tooltip>
      <Button small icon="trash" onClick={onRemove} />
    </div>
  )

}
