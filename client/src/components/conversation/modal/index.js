
import React, { useRef, useState, useEffect } from "react";
import { useDispatch, useSelector } from "react-redux";

import { Input, Modal } from "antd";
import Draggable from "react-draggable";
import cls from "./index.css";

import { findConversationRecord } from "../../../reducers/conversation";
import { ConversationRecord } from "../record";
import { messageOutputSource } from "../util";

const { TextArea } = Input;

const defaultBounds = { left: 0, top: 0, bottom: 0, right: 0 };

const ConversationModal = () => {
  const dispatch = useDispatch();
  const { conversation, currentConversationRecordId } = useSelector(
    (state) => ({
      conversation: state.conversation.conversation,
      currentConversationRecordId:
        state.conversation.currentConversationRecordId,
    })
  );

  // 当前会话（排除根会话）
  const currentConversationRecord = (() => {
    if (
      !currentConversationRecordId ||
      conversation.id === currentConversationRecordId
    ) {
      return null;
    }
    const [record] = findConversationRecord(
      conversation,
      (r) => r.type === "conversation" && r.id === currentConversationRecordId
    );
    return record;
  })();

  const draggleRef = useRef(null);
  const [bounds, setBounds] = useState(defaultBounds);

  const recordsRef = useRef(null);
  const [input, setInput] = useState("");
  const [output, setOutput] = useState("");

  useEffect(() => {
    const ele = recordsRef.current;
    if (!ele) {
      return;
    }
    ele.scrollTop = ele.scrollHeight;
  }, [currentConversationRecord?.data.records.length, output]);

  if (!currentConversationRecord) {
    return null;
  }

  const closeModal = () => {
    dispatch({
      type: "conversation: update current conversation record",
      recordId: null,
    });
  };

  const onStart = (_event, uiData) => {
    const { clientWidth, clientHeight } = window.document.documentElement;
    const targetRect = draggleRef.current?.getBoundingClientRect();
    if (!targetRect) {
      return;
    }
    setBounds({
      left: -targetRect.left + uiData.x,
      right: clientWidth - (targetRect.right - uiData.x),
      top: -targetRect.top + uiData.y,
      bottom: clientHeight - (targetRect.bottom - uiData.y),
    });
  };

  const onChange = (e) => setInput(e.target.value);

  const sendMessage = () => {
    if (input.trim().length === 0 || output.length > 0) {
      return;
    }
    dispatch({
      type: "conversation: add message in conversation record",
      data: {
        recordId: currentConversationRecordId,
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
          type: "conversation: add message in conversation record",
          data: {
            recordId: currentConversationRecordId,
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

  if (!currentConversationRecord) {
    return null;
  }

  const modalRender = (modal) => (
    <Draggable
      bounds={bounds}
      nodeRef={draggleRef}
      onStart={onStart}
      defaultClassName={cls.mask}
    >
      <div ref={draggleRef}>{modal}</div>
    </Draggable>
  );

  const { records } = currentConversationRecord.data;

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
    <Modal
      open
      onCancel={closeModal}
      mask={false}
      maskClosable={false}
      classNames={{ mask: cls.mask, wrapper: cls.mask, body: cls.root }}
      modalRender={modalRender}
      footer={null}
      title={`Conversation(${currentConversationRecordId})`}
    >
      <div ref={recordsRef} className={cls.records}>
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
    </Modal>
  );
};

export default ConversationModal;
