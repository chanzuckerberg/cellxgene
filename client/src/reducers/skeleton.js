const Skeleton = (
  /*
  Skeleton reducer, holds state related to loading/skeleton functionality. 
   */
  state = {
    // data loading flag
    loading: true,
    error: null,

    init: false, // TODO(cc) temp workaround for skeleton PoC first pass, used to skip skeleton before initial load
    skeleton: false,
  },
  action
) => {
  switch (action.type) {
    case "initial data load start":
      return {
        ...state,
        loading: true,
        error: null,
        skeleton: state.init, // Skeleton visible only after initial load (see TODO above)
      };
    case "initial data load complete":
      return {
        ...state,
        loading: false,
        error: null,
        init: true,
        skeleton: false,
      };
    case "initial data load error":
      return {
        ...state,
        error: action.error,
        init: false,
        skeleton: false,
      };
    default:
      return state;
  }
};

export default Skeleton;
