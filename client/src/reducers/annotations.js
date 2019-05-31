// jshint esversion: 6
const Annotations = (
  state = {
    isEditingCategoryName: false,
    isEditingLabelName: false,
    categoryEditable: { category: null },
    labelEditable: { category: null, label: null }
  },
  action
) => {
  switch (action.type) {
    /* CATEGORY */
    case "duplicate annotation category":
      console.log("duplicate category as editable field", action);
      return {
        ...state
      };
    case "add new label to category":
      console.log("add new label to category", action);
      return {
        ...state
        /* Bruce to persist new category name */
      };
    case "activate category edit mode":
      console.log("activate category edit mode", action);
      return {
        ...state,
        isEditingCategoryName: true,
        categoryEditable: {
          category: action.metadataField
        }
      };
    case "category edited":
      console.log("activate category edit mode", action);
      return {
        ...state,
        isEditingCategoryName: true,
        categoryEditable: {
          category: null
        }
        /* Bruce to persist new category name */
      };
    case "delete category":
      console.log("delete category", action);
      return {
        ...state
        /* Bruce to delete category */
      };
    /* LABEL */
    case "activate edit label mode":
      console.log("edit label mode", action);
      return {
        ...state,
        isEditingLabelName: true,
        labelEditable: {
          category: action.metadataField,
          label: action.categoryIndex
        }
      };
    case "cancel edit label mode":
      console.log("cancel edit label mode", action);
      return {
        ...state,
        isEditingLabelName: false,
        labelEditable: { category: null, label: null }
      };
    case "label edited":
      console.log("label edited", action);
      return {
        ...state,
        isEditingLabelName: false,
        labelEditable: { category: null, label: null }
        /* Bruce to persist new label name */
      };
    case "label current cell selection":
      console.log("label cells", action);

      return {
        ...state
        /* add currently selected cells to this label */
      };
    case "delete label":
      console.log("delete label", action);
      return {
        ...state
        /* delete label, add cells to unassigned */
      };
    default:
      return state;
  }
};

export default Annotations;
