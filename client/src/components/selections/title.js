import React from 'react'
import { Button } from '@blueprintjs/core'

export const SelectionGroupTitle = ({ name, onRemove }) => (
  <div style={{ display: 'flex', alignItems: 'center', justifyContent: 'space-between', gap: '1ch' }}>
    <span>{name}</span>
    <Button small icon="trash" onClick={onRemove} />
  </div>
)
