import React from "react";
import cls from "./index.css";
import { RecordTag } from "./tag";

const RoleNameMap = {
  user: "User Question",
  assistant: "Agent Response",
};

export const ConversationMessage = ({ message }) => {
  const { role, content, ctime } = message;

  const name = RoleNameMap[role] || role;

  return (
    <div className={cls.message}>
      <div className={cls.header}>
        <RecordTag tag={name} />
        <span className={cls.ctime}>{new Date(ctime).toLocaleString()}</span>
      </div>
      <span>{content}</span>
    </div>
  );
};
