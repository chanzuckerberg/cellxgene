import { nanoid } from "nanoid";

// type RecordType = "message" | "conversation" | "log";

// interface Message {
//   id: string;
//   role: "user" | "assistant" | "system",
//   content: string;
//   ctime: number;
// }

// interface Log {
//   id: string;
//   type: string;
//   content: string;
//   ctime: number;
// }

// interface Conversation {
//   id: string;
//   records: Record[];
//   ctime: number;
// }

// interface Record {
//   id: string;
//   type: RecordType;
//   data: Message | Conversation | Log;
// }

export const createMessage = ({ role, content }) => ({
  id: `message-${nanoid(6)}`,
  role,
  content,
  ctime: Date.now(),
});

export const createLog = ({ type, content }) => ({
  id: `log-${nanoid(6)}`,
  type,
  content,
  ctime: Date.now(),
});

export const createConversation = () => ({
  id: `conversation-${nanoid(6)}`,
  records: [],
  ctime: Date.now(),
});

const recordCreator = {
  message: createMessage,
  conversation: createConversation,
  log: createLog,
};

export const createRecord = (type, data) => {
  const creator = recordCreator[type];
  if (!creator) {
    throw new Error(`Unknown Record Type: ${type}`);
  }
  const recordData = creator(data);
  return {
    id: `record-${nanoid(6)}`,
    type,
    data: recordData,
  };
};

/**
 * 从root中找到指定conversation
 * @param {*} root 根会话
 * @param {*} conversationId 目标会话ID
 * @returns Conversation | undefined
 */
export const getConversationById = (root, conversationId) => {
  if (root.id === conversationId) {
    return root;
  }
  return root.records.find(
    (r) => r.type === "conversation" && r.data.id === conversationId
  );
};

const initialState = {
  conversation: createConversation(),
};

const ConversationReducer = (state = initialState, action) => {
  switch (action.type) {
    case "conversation: add log": {
      const { type, content } = action;
      const record = createRecord("log", { type, content });
      return {
        ...state,
        conversation: {
          ...state.conversation,
          records: [...state.conversation.records, record],
        },
      };
    }

    case "conversation: add conversation": {
      const record = createRecord("conversation");
      return {
        ...state,
        conversation: {
          ...state.conversation,
          records: [...state.conversation.records, record],
        },
      };
    }

    case "conversation: add message": {
      const { role, content } = action;
      const record = createRecord("message", { role, content });
      return {
        ...state,
        conversation: {
          ...state.conversation,
          records: [...state.conversation.records, record],
        },
      };
    }

    default: {
      return state;
    }
  }
};

export default ConversationReducer;
