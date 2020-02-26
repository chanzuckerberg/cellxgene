import React from "react";
import { connect } from "react-redux";
import {
  Button,
  Menu,
  MenuItem,
  Popover,
  Position,
  PopoverInteractionKind
} from "@blueprintjs/core";

@connect(state => ({
  annotations: state.annotations
}))
class AnnoMenuCategory extends React.Component {
  constructor(props) {
    super(props);
    this.state = {};
  }

  activateAddNewLabelMode = () => {
    const { dispatch, metadataField } = this.props;
    dispatch({
      type: "annotation: activate add new label mode",
      data: metadataField
    });
  };

  activateEditCategoryMode = () => {
    const { dispatch, metadataField } = this.props;

    dispatch({
      type: "annotation: activate category edit mode",
      data: metadataField
    });
  };

  handleDeleteCategory = () => {
    const { dispatch, metadataField } = this.props;
    dispatch({
      type: "annotation: delete category",
      metadataField
    });
  };

  render() {
    const {
      metadataField,
      annotations,
      isUserAnno,
      createText,
      editText,
      deleteText
    } = this.props;

    return (
      <>
        {isUserAnno ? (
          <Popover
            interactionKind={PopoverInteractionKind.HOVER}
            boundary="window"
            position={Position.RIGHT_TOP}
            content={
              <Menu>
                <MenuItem
                  icon="tag"
                  data-testclass="handleAddNewLabelToCategory"
                  data-testid={`${metadataField}:add-new-label-to-category`}
                  onClick={this.activateAddNewLabelMode}
                  text={createText}
                />
                <MenuItem
                  icon="edit"
                  disabled={annotations.isEditingCategoryName}
                  data-testclass="activateEditCategoryMode"
                  data-testid={`${metadataField}:edit-category-mode`}
                  onClick={this.activateEditCategoryMode}
                  text={editText}
                />
                <MenuItem
                  icon="delete"
                  intent="danger"
                  data-testclass="handleDeleteCategory"
                  data-testid={`${metadataField}:delete-category`}
                  onClick={this.handleDeleteCategory}
                  text={deleteText}
                />
              </Menu>
            }
          >
            <Button
              style={{ marginLeft: 0 }}
              data-testclass="seeActions"
              data-testid={`${metadataField}:see-actions`}
              icon="more"
              minimal
            />
          </Popover>
        ) : null}
      </>
    );
  }
}

export default AnnoMenuCategory;
