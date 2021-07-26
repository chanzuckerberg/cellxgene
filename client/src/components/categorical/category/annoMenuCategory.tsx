import React from "react";
import { connect } from "react-redux";
import {
  Button,
  Menu,
  MenuItem,
  Popover,
  Position,
  Tooltip,
  Icon,
  PopoverInteractionKind,
  Intent,
} from "@blueprintjs/core";
import { IconNames } from "@blueprintjs/icons";

import * as globals from "../../../globals";
import actions from "../../../actions";

type State = any;

// @ts-expect-error ts-migrate(1238) FIXME: Unable to resolve signature of class decorator whe... Remove this comment to see the full error message
@connect((state) => ({
  annotations: (state as any).annotations,
}))
class AnnoMenuCategory extends React.PureComponent<{}, State> {
  constructor(props: {}) {
    super(props);
    this.state = {};
  }

  activateAddNewLabelMode = () => {
    // @ts-expect-error ts-migrate(2339) FIXME: Property 'dispatch' does not exist on type 'Readon... Remove this comment to see the full error message
    const { dispatch, metadataField } = this.props;
    dispatch({
      type: "annotation: activate add new label mode",
      data: metadataField,
    });
  };

  activateEditCategoryMode = () => {
    // @ts-expect-error ts-migrate(2339) FIXME: Property 'dispatch' does not exist on type 'Readon... Remove this comment to see the full error message
    const { dispatch, metadataField } = this.props;
    dispatch({
      type: "annotation: activate category edit mode",
      data: metadataField,
    });
  };

  handleDeleteCategory = () => {
    // @ts-expect-error ts-migrate(2339) FIXME: Property 'dispatch' does not exist on type 'Readon... Remove this comment to see the full error message
    const { dispatch, metadataField } = this.props;
    dispatch(actions.annotationDeleteCategoryAction(metadataField));
  };

  render() {
    const {
      // @ts-expect-error ts-migrate(2339) FIXME: Property 'metadataField' does not exist on type 'R... Remove this comment to see the full error message
      metadataField,
      // @ts-expect-error ts-migrate(2339) FIXME: Property 'annotations' does not exist on type 'Rea... Remove this comment to see the full error message
      annotations,
      // @ts-expect-error ts-migrate(2339) FIXME: Property 'isUserAnno' does not exist on type 'Read... Remove this comment to see the full error message
      isUserAnno,
      // @ts-expect-error ts-migrate(2339) FIXME: Property 'createText' does not exist on type 'Read... Remove this comment to see the full error message
      createText,
      // @ts-expect-error ts-migrate(2339) FIXME: Property 'editText' does not exist on type 'Readon... Remove this comment to see the full error message
      editText,
      // @ts-expect-error ts-migrate(2339) FIXME: Property 'deleteText' does not exist on type 'Read... Remove this comment to see the full error message
      deleteText,
    } = this.props;
    return (
      <>
        {isUserAnno ? (
          <>
            <Tooltip
              content={createText}
              position="bottom"
              hoverOpenDelay={globals.tooltipHoverOpenDelay}
            >
              <Button
                style={{ marginLeft: 0, marginRight: 2 }}
                data-testclass="handleAddNewLabelToCategory"
                data-testid={`${metadataField}:add-new-label-to-category`}
                icon={<Icon icon="plus" iconSize={10} />}
                onClick={this.activateAddNewLabelMode}
                small
                minimal
              />
            </Tooltip>
            <Popover
              interactionKind={PopoverInteractionKind.HOVER}
              boundary="window"
              position={Position.RIGHT_TOP}
              content={
                <Menu>
                  <MenuItem
                    icon="edit"
                    disabled={annotations.isEditingCategoryName}
                    data-testclass="activateEditCategoryMode"
                    data-testid={`${metadataField}:edit-category-mode`}
                    onClick={this.activateEditCategoryMode}
                    text={editText}
                  />
                  <MenuItem
                    icon={IconNames.TRASH}
                    intent={Intent.DANGER}
                    data-testclass="handleDeleteCategory"
                    data-testid={`${metadataField}:delete-category`}
                    onClick={this.handleDeleteCategory}
                    text={deleteText}
                  />
                </Menu>
              }
            >
              <Button
                style={{ marginLeft: 0, marginRight: 5 }}
                data-testclass="seeActions"
                data-testid={`${metadataField}:see-actions`}
                icon={<Icon icon="more" iconSize={10} />}
                small
                minimal
              />
            </Popover>
          </>
        ) : null}
      </>
    );
  }
}

export default AnnoMenuCategory;
