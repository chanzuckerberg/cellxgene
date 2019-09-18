/*
Reducers for annotation UI-state.
*/
const Annotations = (
  state = {
    isEditingCategoryName: false,
    isEditingLabelName: false,
    categoryEditable: false,
    labelEditable: { category: null, label: null }
  },
  action
) => {
  switch (action.type) {
    /* CATEGORY */
    case "annotation: activate add new label mode":
      console.log(action.type, action);
      return {
        ...state,
        isAddingNewLabel: true,
        categoryAddingNewLabel: action.data
      };
    case "annotation: disable add new label mode":
      return {
        ...state,
        isAddingNewLabel: false,
        categoryAddingNewLabel: null
      };
    case "annotation: add new label to category":
      console.log(action.type, action);
      return {
        ...state,
        isAddingNewLabel: false,
        categoryAddingNewLabel: null
      };
    case "annotation: activate category edit mode":
      console.log(action.type, action);
      return {
        ...state,
        isEditingCategoryName: true,
        categoryEditable: action.data
      };
    case "annotation: disable category edit mode":
      console.log(action.type, action);
      return {
        ...state,
        isEditingCategoryName: false,
        categoryEditable: null
      };
    case "annotation: category edited":
      console.log(action.type, action);
      return {
        ...state,
        isEditingCategoryName: true,
        categoryEditable: null
      };

    /* LABEL */
    case "annotation: activate edit label mode":
      console.log(action.type, action);
      return {
        ...state,
        isEditingLabelName: true,
        labelEditable: {
          category: action.metadataField,
          label: action.categoryIndex
        }
      };
    case "annotation: cancel edit label mode":
      console.log(action.type, action);
      return {
        ...state,
        isEditingLabelName: false,
        labelEditable: { category: null, label: null }
      };
    case "annotation: label edited":
      console.log(action.type, action);
      return {
        ...state,
        isEditingLabelName: false,
        labelEditable: { category: null, label: null }
        /* Bruce to persist new label name */
      };
    default:
      return state;
  }
};

export default Annotations;
