/*
Reducers for annotation UI-state.
*/
const Annotations = (
  state = {
    /* 
    Annotations collection name - which will be used to save the named set of annotations
    in some persistent store (database, file system, etc).

    Backend may expect this to be a legal file name, which is typically alpha-numeric, plus [_-,.].
    Keep it simple.

    If `annotationDataCollectionNameIsReadOnly` is true, you may NOT change the data collection name.
    If it is true, you may change it and it will be used at the time the annotations are written to the
    back-end.
    */
    dataCollectionNameIsReadOnly: true,
    dataCollectionName: null,

    /* UI component state */
    isEditingCategoryName: false,
    isEditingLabelName: false,
    categoryBeingEdited: null,
    categoryAddingNewLabel: null,
    labelEditable: { category: null, label: null }
  },
  action
) => {
  switch (action.type) {
    case "configuration load complete": {
      const dataCollectionName =
        action.config.parameters?.["annotations-data-collection-name"] ?? null;
      const dataCollectionNameIsReadOnly =
        action.config.parameters?.[
          "annotations-data-collection-name-is-read-only"
        ] ?? false;
      return {
        ...state,
        dataCollectionNameIsReadOnly,
        dataCollectionName
      };
    }

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
