/**
 * // 消息类型
 * interface Message {
 *   // 消息ID
 *   id: string;
 *   // 角色
 *   role: UserRole;
 *   // 消息内容
 *   content: string;
 *   // 消息时间（时间戳）
 *   timestamp: number;
 * }
 *
 * // 选区类型
 * interface Selection {
 *   // 嵌入名称
 *   emb: string;
 *   // 选区名称
 *   name: string;
 *   // 多边形顶点坐标列表
 *   polygon: vec2[]
 *   // 选中的细胞Index
 *   indexes: number[]
 * };
 *
 * // 会话上下文
 * interface SessionContext {
 *  // 选区列表
 *  selections: Selection[];
 * }
 *
 * // 会话类型
 * interface Session {
 *  // 会话上下文
 *  context: SessionContext;
 *  // 消息列表
 *  messages: Message[];
 * }
 */

const initialState = {
  session: {
    context: {
      selections: [],
    },
    messages: [],
  },
};

const ChatReducer = (state = initialState, action) => {
  switch (action.type) {
    case "chat: new message": {
      const { message } = action;
      return {
        ...state,
        session: {
          ...state.session,
          messages: [...state.session.messages, message],
        },
      };
    }
    default: {
      return state;
    }
  }
};

export default ChatReducer;
