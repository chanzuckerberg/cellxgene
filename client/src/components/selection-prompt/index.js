import React, { useEffect, useState } from "react";
import { Dialog, Classes, TextArea } from "@blueprintjs/core";

const SelectionPrompt = ({ visible = false, onClose, onSubmit }) => {
  const [q, setQ] = useState("");

  const submit = async () => {
    onSubmit(q);
  };

  useEffect(() => {
    if (visible) {
      setQ("");
    }
  }, [visible]);

  return (
    <Dialog isOpen={visible} onClose={onClose}>
      <div className={Classes.DIALOG_BODY}>
        <TextArea
          autoFocus
          fill
          style={{ resize: "none" }}
          growVertically
          placeholder="Ask a question related to this selection，Press Enter to send, Shift+Enter for newline"
          value={q}
          onChange={(e) => setQ(e.target.value)}
          onKeyDown={async (e) => {
            if (e.key === "Enter" && !e.shiftKey) {
              e.stopPropagation();
              e.preventDefault();
              await submit();
            }
          }}
        />
      </div>
    </Dialog>
  );
};

export default SelectionPrompt;
