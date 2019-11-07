/*
Reducers for annotation UI-state.
*/
const Annotations = (
  state = {
    isEditingCategoryName: false,
    isEditingLabelName: false,
    categoryBeingEdited: null,
    categoryAddingNewLabel: null,
    labelEditable: { category: null, label: null }
  },
  action
) => {
  switch (action.type) {
    /* CATEGORY */
    case "annotation: activate add new label mode":
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
      return {
        ...state,
        isAddingNewLabel: false,
        categoryAddingNewLabel: null
      };
    case "annotation: activate category edit mode":
      return {
        ...state,
        isEditingCategoryName: true,
        categoryBeingEdited: action.data
      };
    case "annotation: disable category edit mode":
      return {
        ...state,
        isEditingCategoryName: false,
        categoryBeingEdited: null
      };
    case "annotation: category edited":
      return {
        ...state,
        isEditingCategoryName: false,
        categoryBeingEdited: null
      };

    /* LABEL */
    case "annotation: activate edit label mode":
      return {
        ...state,
        isEditingLabelName: true,
        labelEditable: {
          category: action.metadataField,
          label: action.categoryIndex
        }
      };
    case "annotation: cancel edit label mode":
      return {
        ...state,
        isEditingLabelName: false,
        labelEditable: { category: null, label: null }
      };
    case "annotation: label edited":
      return {
        ...state,
        isEditingLabelName: false,
        labelEditable: { category: null, label: null }
      };
    default:
      return state;
  }
};

export default Annotations;
