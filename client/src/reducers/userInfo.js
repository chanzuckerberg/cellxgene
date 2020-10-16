const UserInfo = (state = {}, action) => {
  switch (action.type) {
    case "initial data load start":
      return {
        ...state,
        loading: true,
        error: null,
      };
    case "userInfo load complete":
      return {
        ...state,
        loading: false,
        error: null,
        ...action.userInfo,
      };
    case "initial data load error":
      return {
        ...state,
        error: action.error,
      };
    default:
      return state;
  }
};

export default UserInfo;
