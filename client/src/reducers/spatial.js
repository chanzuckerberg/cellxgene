const Spatial = (
  state = {
    loading: null,
    error: null,
    metadata: null,
  },
  action
) => {
  switch (action.type) {
    case "request spatial metadata started":
      return {
        ...state,
        loading: true,
        error: null,
      };
    case "request spatial metadata success":
      return {
        ...state,
        error: null,
        loading: false,
        metadata: action,
      };
    case "request spatial metadata error":
      return {
        ...state,
        loading: false,
        error: action.data,
      };
    default:
      return state;
  }
};

export default Spatial;
