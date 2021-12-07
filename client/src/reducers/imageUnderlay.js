const imageUnderlay = (state = { isActive: false }, action) => {
  switch (action.type) {
    case "toggle image underlay":
      return {
        ...state,
        isActive: !state.isActive,
      };

    default:
      return state;
  }
};

export default imageUnderlay;
