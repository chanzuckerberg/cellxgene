import React, { useEffect, useRef, useState } from "react";
import { useDispatch, useSelector } from "react-redux";
import { Button } from "@blueprintjs/core";
import { Input } from "antd";
import { PushpinOutlined } from "@ant-design/icons";
import cls from "./index.css";

import { messageOutputSource } from "./util";
import { ConversationRecord } from "./record";
import ConversationModal from "./modal";

const { TextArea } = Input;

const Conversation = () => {
  const dispatch = useDispatch();
  const { conversation, pinned } = useSelector((state) => ({
    conversation: state.conversation.conversation,
    pinned: state.pin.right,
  }));

  const ref = useRef(null);
  const [input, setInput] = useState("");
  const [output, setOutput] = useState("");

  const archive = () => {
    const text = "test message";
    const blob = new Blob([text], { type: "text/plain" });
    const url = URL.createObjectURL(blob);
    const a = document.createElement("a");
    a.href = url;
    a.download = "chat_archive.txt";
    a.click();
    URL.revokeObjectURL(url);
  };

  const onChange = (e) => setInput(e.target.value);

  const sendMessage = () => {
    if (input.trim().length === 0 || output.length > 0) {
      return;
    }
    dispatch({
      type: "conversation: add message in root",
      data: {
        role: "user",
        content: input,
      },
    });
    setInput("");

    let chunks = "";
    messageOutputSource.subscribe({
      next: (chunk) => {
        chunks += chunk;
        setOutput(chunks);
      },
      complete: () => {
        dispatch({
          type: "conversation: add message in root",
          data: {
            role: "assistant",
            content: chunks,
          },
        });
        setOutput("");
      },
    });
  };

  const onKeyDown = (e) => {
    if (e.key === "Enter" && !e.shiftKey) {
      e.stopPropagation();
      e.preventDefault();
      sendMessage();
    }
  };

  const unpin = () => {
    dispatch({
      type: "pin: update",
      loc: "right",
      pinned: false,
    });
  };

  useEffect(() => {
    const ele = ref.current;
    if (!ele) {
      return;
    }
    ele.scrollTop = ele.scrollHeight;
  }, [conversation.records.length, output]);

  const { records } = conversation;

  const tmpRecord = {
    id: "record-tmp",
    type: "message",
    data: {
      role: "assistant",
      content: output,
      ctime: Date.now(),
    },
  };

  return (
    <>
      <div className={cls.root}>
        <div className={cls.header}>
          {pinned && <Button icon={<PushpinOutlined />} onClick={unpin} />}
          <span>Conversation</span>
          <Button
            icon="archive"
            disabled={records.length === 0 || output.length > 0}
            onClick={archive}
          />
        </div>

        <div ref={ref} className={cls.records}>
          {records.map((r) => (
            <ConversationRecord key={r.id} record={r} />
          ))}
          {tmpRecord.data.content && (
            <ConversationRecord key={tmpRecord.id} record={tmpRecord} />
          )}
        </div>

        <TextArea
          value={input}
          onChange={onChange}
          onKeyDown={onKeyDown}
          placeholder="Press Enter to send, Shift+Enter for newline"
          autoSize={{ minRows: 3, maxRows: 6 }}
        />
      </div>

      <ConversationModal />
    </>
  );
};

export default Conversation;
