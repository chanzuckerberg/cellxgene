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
//   selectionIds: number[];
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
  selectionIds: [],
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
 * 从root records中找到指定record
 * @param {*} root 根会话
 * @param {*} finder 检索函数, `(record, index) => boolean`
 * @returns [record, recordIndex]
 */
export const findConversationRecord = (root, finder) => {
  const recordIndex = root.records.findIndex(finder);
  if (recordIndex === -1) {
    return [null, recordIndex];
  }
  return [root.records[recordIndex], recordIndex];
};

const initialConversation = (() => {
  const c = createConversation();
  // c.records.push(createRecord("message", { role: 'user', content: 'hi' }))
  // c.records.push(createRecord("message", { role: 'assistant', content: 'hello, what can i do?' }))
  // c.records.push(createRecord("log", { type: 'Selection Create', content: '123 cells' }))
  // c.records.push(createRecord("conversation"))
  // c.records.push(createRecord("log", { type: 'Selection Remove', content: '123 cells' }))
  return c;
})();

const initialState = {
  conversation: initialConversation,
  currentConversationRecordId: null,
};

const ConversationReducer = (state = initialState, action) => {
  switch (action.type) {
    case "conversation: add log": {
      const { data } = action;
      const { type, content } = data;
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
        currentConversationRecordId: record.id,
        conversation: {
          ...state.conversation,
          records: [...state.conversation.records, record],
        },
      };
    }

    case "conversation: add message in root": {
      // 直接在根会话添加消息
      const { data } = action;
      const { role, content } = data;

      const record = createRecord("message", { role, content });

      return {
        ...state,
        conversation: {
          ...state.conversation,
          records: [...state.conversation.records, record],
        },
      };
    }

    case "conversation: add message in conversation record": {
      // 在会话记录中添加消息
      const { data } = action;
      const { recordId, role, content } = data;

      const record = createRecord("message", { role, content });

      // 找到目标会话记录
      const [conversationRecord, recordIndex] = findConversationRecord(
        state.conversation,
        (r) => r.type === "conversation" && r.id === recordId
      );
      if (recordIndex < 0) {
        return state;
      }
      return {
        ...state,
        conversation: {
          ...state.conversation,
          records: [
            ...state.conversation.records.slice(0, recordIndex),
            {
              ...conversationRecord,
              data: {
                ...conversationRecord.data,
                records: [...conversationRecord.data.records, record],
              },
            },
            ...state.conversation.records.slice(recordIndex + 1),
          ],
        },
      };
    }

    case "conversation: update current conversation record": {
      const { recordId } = action;

      return {
        ...state,
        currentConversationRecordId: recordId,
      };
    }

    default: {
      return state;
    }
  }
};

export default ConversationReducer;
