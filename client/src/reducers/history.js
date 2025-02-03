export const ActionType = {
  // Selection
  SELECTION_CREATE: "SELECTION_CREATE",
  SELECTION_SET_CREATE: "SELECTION_SET_CREATE",
  SELECTION_UPDATE: "SELECTION_UPDATE",
  SELECTION_DELETE: "SELECTION_DELETE",
  // Chat
  CHAT_SEND: "CHAT_SEND",
  CHAT_RECEIVE: "CHAT_RECEIVE",
};

// interface Action {
//   type: ActionType[keyof ActionType],
//   payload: Object
// }

const initialState = {
  actions: [],
};

const HistoryReducer = (state = initialState, action) => {
  switch (action.type) {
    default: {
      return state;
    }
  }
};

export default HistoryReducer;
