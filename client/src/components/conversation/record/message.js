import React from "react";
import cls from "./index.css";


export const ConversationMessage = ({ message }) => {
  const { role, content, ctime } = message;

  return (
    <div className={cls.message}>
      <div className={cls.header}>
        <span className={cls.role}>{role}</span>
        <span className={cls.ctime}>{new Date(ctime).toLocaleString()}</span>
      </div>
      <span className={cls.content}>{content}</span>
    </div>
  );
};
