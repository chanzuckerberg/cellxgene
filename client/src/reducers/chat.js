/**
 * // 消息类型
 * interface Message {
 *   // 消息ID
 *   id: number;
 *   // 用户ID
 *   userId:  number;
 *   // 角色
 *   role: UserRole;
 *   // 消息类型
 *   content: string;
 *   // 消息时间（时间戳）
 *   timestamp: number;
 * }
 *
 * // 选区类型
 * interface Selection {
 *   // 选区ID
 *   id: string;
 *   // 多边形顶点坐标列表
 *   polygon: vec2[]
 *   // 数据集ID
 *   datasetId: string;
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
 *  // 会话ID
 *  id: number;
 *  // 会话上下文
 *  context: SessionContext;
 *  // 消息列表
 *  messages: Message[];
 * }
 *
 * // 临时聊天会话
 * interface ChatSessionTemp {
 *  // 输入框是否可见
 *  visible: string;
 *  // 内容
 *  content: string;
 * }
 *
 * // 聊天状态
 * interface ChatState {
 *  // 会话列表
 *  sessions: Session[];
 *  // 临时状态
 *  tmp: ChatSessionTemp;
 * }
 */

const DefaultState = {
  sessions: [],
  tmp: {
    visible: false,
    content: "",
  },
};

const CACHE_KEY = "chat_state";

// 初始化时从 localStorage 读取状态
const loadState = () => {
  try {
    const serializedState = localStorage.getItem(CACHE_KEY);
    if (serializedState === null) {
      return DefaultState;
    }
    return JSON.parse(serializedState);
  } catch (err) {
    console.error("Failed to load state from localStorage", err);
    return DefaultState;
  }
};

// 保存状态到 localStorage
// const saveState = (state) => {
//   try {
//     const serializedState = JSON.stringify(state);
//     localStorage.setItem(CACHE_KEY, serializedState);
//   } catch (err) {
//     console.error("Failed to save state to localStorage", err);
//   }
// };
const saveState = () => {
  // 暂时不保存
  
};

const initialState = loadState();

const ChatReducer = (state = initialState, action) => {
  switch (action.type) {
    // case "chat: load state": {
    //   return loadState()
    // }

    // case 'chat: save state': {
    //   saveState()
    //   return state
    // }

    case "chat: open temp session": {
      return {
        ...state,
        tmp: {
          ...state.tmp,
          visible: true,
          content: "",
        },
      };
    }

    case "chat: close temp session": {
      return {
        ...state,
        tmp: {
          ...state.tmp,
          visible: false,
        },
      };
    }

    case "chat: create session": {
      console.log(action);
      const { prompt, polygon, indexes } = action;
      const session = {
        id: 0,
        context: {
          selections: [
            {
              id: "",
              // 多边形顶点坐标列表
              polygon,
              // 数据集ID
              datasetId: "",
              // 选中的细胞Index
              indexes,
            },
          ],
        },
        messages: [
          {
            id: "",
            userId: 0,
            role: "user",
            content: prompt,
            timestamp: Date.now(),
          },
        ],
      };
      const newState = {
        ...state,
        tmp: {
          visible: false,
          content: "",
        },
        sessions: [...state.sessions, session],
      };

      saveState(newState);

      return newState;
    }

    default: {
      return state;
    }
  }
};

export default ChatReducer;
