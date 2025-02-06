
import { Tag } from "antd";
import React from "react";
import cls from "./index.css";

const TagColorMap = {
  "User Question": "#2beca5",
  "Agent Response": "#2bc9ec",
  Conversation: "#ec322b",
  "Selection Create": "#d377ec",
  "Selection Remove": "#ec77a8",
};

export const RecordTag = ({ tag }) => {
  const color = TagColorMap[tag];

  return (
    <Tag className={cls.tag} color={color}>
      {tag}
    </Tag>
  );
};
