import React from "react";
import cls from "./index.css";


import { RecordTag } from "./tag";

export const ConversationLog = ({ log }) => {
  const { type, ctime } = log;

  return (
    <div className={cls.log}>
      <RecordTag tag={type} />
      <span className={cls.ctime}>{new Date(ctime).toLocaleString()}</span>
    </div>
  );
};
