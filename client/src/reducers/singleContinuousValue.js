const initialState = {
  singleContinuousValues: new Map(),
};
const singleContinuousValue = (state = initialState, action) => {
  switch (action.type) {
    case "add single continuous value":
      state.singleContinuousValues.set(action.field, action.value);
      return state;
    default:
      return state;
  }
};

export default singleContinuousValue;
