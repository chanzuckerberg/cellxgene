/*
Reducers for annotation UI-state.
*/
const Annotations = (
  state = {
    /*
    Annotations collection name - which will be used to save the named set of annotations
    in some persistent store (database, file system, etc).

    The backend may expect this to be a legal file name, which is typically alpha-numeric, plus [_-,.].
    Keep it simple or the server may return an error.

    If `dataCollectionNameIsReadOnly` is true, you may NOT change the data collection name.
    If false, you may change `dataCollectionName` and it will be used at the time the annotations are
    written to the back-end.

    */
    dataCollectionNameIsReadOnly: true,
    dataCollectionName: null,

    /*
    Annotations UI component state
    */
    isEditingCategoryName: false,
    isEditingLabelName: false,
    categoryBeingEdited: null,
    categoryAddingNewLabel: null,
    labelEditable: { category: null, label: null },
    promptForFilename: true,
  },
  action: any
) => {
  switch (action.type) {
    case "configuration load complete": {
      const dataCollectionName =
        action.config.parameters?.["annotations-data-collection-name"] ?? null;
      const dataCollectionNameIsReadOnly =
        (action.config.parameters?.[
          "annotations-data-collection-is-read-only"
        ] ??
          false) &&
        (action.config.parameters?.annotations_genesets_name_is_read_only ??
          true);

      const promptForFilename =
        action.config.parameters?.user_annotation_collection_name_enabled;
      return {
        ...state,
        dataCollectionNameIsReadOnly,
        dataCollectionName,
        promptForFilename,
      };
    }

    case "set annotations collection name": {
      if (state.dataCollectionNameIsReadOnly) {
        throw new Error("data collection name is read only");
      }
      return {
        ...state,
        dataCollectionName: action.data,
      };
    }

    /* CATEGORY */
    case "annotation: activate add new label mode": {
      return {
        ...state,
        isAddingNewLabel: true,
        categoryAddingNewLabel: action.data,
      };
    }

    case "annotation: disable add new label mode": {
      return {
        ...state,
        isAddingNewLabel: false,
        categoryAddingNewLabel: null,
      };
    }

    case "annotation: activate category edit mode": {
      return {
        ...state,
        isEditingCategoryName: true,
        categoryBeingEdited: action.data,
      };
    }

    case "annotation: disable category edit mode": {
      return {
        ...state,
        isEditingCategoryName: false,
        categoryBeingEdited: null,
      };
    }

    /* LABEL */
    case "annotation: activate edit label mode": {
      return {
        ...state,
        isEditingLabelName: true,
        labelEditable: {
          category: action.metadataField,
          label: action.categoryIndex,
        },
      };
    }

    case "annotation: cancel edit label mode": {
      return {
        ...state,
        isEditingLabelName: false,
        labelEditable: { category: null, label: null },
      };
    }

    default:
      return state;
  }
};

export default Annotations;
