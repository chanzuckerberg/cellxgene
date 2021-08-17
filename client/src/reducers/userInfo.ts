import { Action } from "redux";

export interface UserInfoAction extends Action<string>, User {
  userInfo: UserInfoPayload;
  error: string;
}

export interface UserInfoPayload {
  is_authenticated: boolean;
  username: string;
  user_id: string;
  email: string;
  picture: string;
}

export interface UserInfoState extends UserInfoPayload {
  loading: boolean;
  error: string | null;
}

const UserInfo = (
  state: UserInfoState,
  action: UserInfoAction
): UserInfoState => {
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
