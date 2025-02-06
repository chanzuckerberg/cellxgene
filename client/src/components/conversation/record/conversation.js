import React from "react";
import { Button } from "antd";
import cls from "./index.css";

import { RecordTag } from "./tag";

export const ConversationConversation = ({ conversation }) => {
  const { ctime } = conversation;

  return (
    <div className={cls.message}>
      <div className={cls.header}>
        <RecordTag tag="Conversation" />
        <span className={cls.ctime}>{new Date(ctime).toLocaleString()}</span>
      </div>
      <Button onClick={() => console.log(conversation)}>Reopen</Button>
    </div>
  );
};
