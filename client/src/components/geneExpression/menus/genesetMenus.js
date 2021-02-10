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
} from "@blueprintjs/core";

import * as globals from "../../../globals";
import AddGeneToGenesetDialogue from "./addGeneToGenesetDialogue";

@connect((state) => {
  return {
    genesetsUI: state.genesetsUI,
  };
})
class GenesetMenus extends React.PureComponent {
  constructor(props) {
    super(props);
    this.state = {};
  }

  activateAddGeneToGenesetMode = () => {
    const { dispatch, geneset } = this.props;
    dispatch({
      type: "geneset: activate add new genes mode",
      geneset,
    });
  };

  activateEditGenesetNameMode = () => {
    const { dispatch, geneset } = this.props;

    dispatch({
      type: "geneset: activate geneset edit mode",
      data: geneset,
    });
  };

  handleDeleteCategory = () => {
    const { dispatch, geneset } = this.props;
    dispatch({ type: "geneset: delete", name: geneset });
  };

  render() {
    const { geneset, genesetsEditable, createText } = this.props;

    return (
      <>
        {genesetsEditable ? (
          <>
            <Tooltip
              content={createText}
              position={Position.BOTTOM}
              hoverOpenDelay={globals.tooltipHoverOpenDelay}
            >
              <Button
                style={{ marginLeft: 0, marginRight: 2 }}
                data-testclass="handleAddNewLabelToCategory"
                data-testid={`${geneset}:add-new-label-to-category`}
                icon={<Icon icon="plus" iconSize={10} />}
                onClick={this.activateAddGeneToGenesetMode}
                small
                minimal
              />
            </Tooltip>
            <AddGeneToGenesetDialogue geneset={geneset} />
            <Popover
              interactionKind={PopoverInteractionKind.HOVER}
              boundary="window"
              position={Position.BOTTOM}
              content={
                <Menu>
                  <MenuItem
                    icon="edit"
                    data-testclass="activateEditGenesetNameMode"
                    data-testid={`${geneset}:edit-genesetName-mode`}
                    onClick={this.activateEditGenesetNameMode}
                    text="Edit gene set name"
                  />
                  <MenuItem
                    icon="trash"
                    intent="danger"
                    data-testclass="handleDeleteCategory"
                    data-testid={`${geneset}:delete-category`}
                    onClick={this.handleDeleteCategory}
                    text="Delete this gene set (destructive, will remove set and collection of genes)"
                  />
                </Menu>
              }
            >
              <Button
                style={{ marginLeft: 0, marginRight: 5 }}
                data-testclass="seeActions"
                data-testid={`${geneset}:see-actions`}
                icon={<Icon icon="more" iconSize={10} />}
                small
                minimal
              />
            </Popover>
            <Button
              style={{ marginLeft: 0 }}
              data-testclass="colorby-entire-geneset"
              data-testid={`${geneset}:colorby-entire-geneset`}
              icon={<Icon icon="tint" iconSize={16} />}
            />
          </>
        ) : null}
      </>
    );
  }
}

export default GenesetMenus;
