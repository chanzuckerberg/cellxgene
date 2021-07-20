const Config = (
  state = {
    displayNames: null,
    features: null,
    parameters: null,
  },
  action: any
) => {
  switch (action.type) {
    case "initial data load start":
      return {
        ...state,
        loading: true,
        error: null,
      };
    case "configuration load complete":
      return {
        ...state,
        loading: false,
        error: null,
        ...action.config,
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

export default Config;
