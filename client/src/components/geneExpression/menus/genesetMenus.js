import React from "react";
import { connect } from "react-redux";
import { Tooltip2 } from "@blueprintjs/popover2";

import {
  Button,
  AnchorButton,
  Menu,
  MenuItem,
  Popover,
  Position,
  Icon,
  PopoverInteractionKind,
} from "@blueprintjs/core";

import * as globals from "../../../globals";
import actions from "../../../actions";
import AddGeneToGenesetDialogue from "./addGeneToGenesetDialogue";

@connect((state) => {
  return {
    genesetsUI: state.genesetsUI,
    colorAccessor: state.colors.colorAccessor,
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
      type: "geneset: activate rename geneset mode",
      data: geneset,
    });
  };

  handleColorByEntireGeneset = () => {
    const { dispatch, geneset } = this.props;

    dispatch({
      type: "color by geneset mean expression",
      geneset,
    });
  };

  handleDeleteCategory = () => {
    const { dispatch, geneset } = this.props;
    dispatch(actions.genesetDelete(geneset));
  };

  render() {
    const {
      geneset,
      genesetsEditable,
      createText,
      colorAccessor,
      isOpen,
      toggleSummaryHisto,
    } = this.props;

    const isColorBy = geneset === colorAccessor;
    const genesetClosed = !isOpen;
    const showingAllGenes = !toggleSummaryHisto;

    const notShowingSummary = genesetClosed || showingAllGenes;

    return (
      <>
        {genesetsEditable && (
          <>
            <Tooltip2
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
            </Tooltip2>
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
                    text="Edit gene set name and description"
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
            <Tooltip2
              content={`Color by gene set ${geneset} mean (gene set must be open, and toggled to mean expression histogram)`}
              position={Position.BOTTOM}
              hoverOpenDelay={globals.tooltipHoverOpenDelay}
            >
              <AnchorButton
                disabled={notShowingSummary && !isColorBy}
                active={isColorBy}
                intent={isColorBy ? "primary" : "none"}
                style={{ marginLeft: 0 }}
                onClick={this.handleColorByEntireGeneset}
                data-testclass="colorby-entire-geneset"
                data-testid={`${geneset}:colorby-entire-geneset`}
                icon={<Icon icon="tint" iconSize={16} />}
              />
            </Tooltip2>
          </>
        )}
      </>
    );
  }
}

export default GenesetMenus;
