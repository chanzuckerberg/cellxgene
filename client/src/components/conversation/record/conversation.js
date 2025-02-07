import React from "react";
import { Button } from "antd";
import { useDispatch } from "react-redux";
import cls from "./index.css";

import { RecordTag } from "./tag";

export const ConversationConversation = ({ recordId, conversation }) => {
  const { ctime } = conversation;

  const dispatch = useDispatch();

  const openConversation = () => {
    dispatch({
      type: "conversation: update current conversation record",
      recordId,
    });
  };

  return (
    <div className={cls.message}>
      <div className={cls.header}>
        <RecordTag tag="Conversation" />
        <span className={cls.ctime}>{new Date(ctime).toLocaleString()}</span>
      </div>
      <Button onClick={openConversation}>Reopen</Button>
    </div>
  );
};
