import React, { useEffect, useState } from "react";
import { Dialog, Classes } from "@blueprintjs/core";
import { Input } from "antd";
import { Button } from '@blueprintjs/core'

export const SelectionCreator = ({ visible, onClose, onCreate }) => {
  const [name, setName] = useState("");

  const submit = () => {
    onCreate({ name });
    onClose();
  };

  useEffect(() => {
    if (visible) {
      setName("");
    }
  }, [visible]);

  return (
    <Dialog isOpen={visible} title="Create Region Set" onClose={onClose}>
      <form
        style={{ display: "flex", alignItems: "center", gap: "2ch" }}
        className={Classes.DIALOG_BODY}
        onSubmit={(e) => {
          e.preventDefault();
          submit();
        }}
      >
        <Input
          style={{ flex: 1 }}
          autoFocus
          placeholder="Region Set Name"
          value={name}
          onChange={(e) => setName(e.target.value)}
        />

        <Button type="submit" intent="primary">
          Create
        </Button>
      </form>
    </Dialog>
  );
};
