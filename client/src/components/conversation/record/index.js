import React from "react";

import { ConversationMessage } from "./message";

export const ConversationRecord = ({ record }) => {
  const { type, data } = record;

  if (type === "message") {
    return <ConversationMessage message={data} />;
  }

  console.warn(
    `unknown conversation record type: ${type}, data: ${JSON.stringify(data)}`
  );
  return null;
};
