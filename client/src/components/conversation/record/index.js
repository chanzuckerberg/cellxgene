import React from "react";

import { ConversationMessage } from "./message";
import { ConversationLog } from "./log";
import { ConversationConversation } from "./conversation";

export const ConversationRecord = ({ record }) => {
  const { type, data } = record;

  if (type === "message") {
    return <ConversationMessage message={data} />;
  }

  if (type === "log") {
    return <ConversationLog log={data} />;
  }

  if (type === "conversation") {
    return <ConversationConversation conversation={data} />;
  }

  console.warn(
    `unknown conversation record type: ${type}, data: ${JSON.stringify(data)}`
  );
  return null;
};
