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

type State = any;

// @ts-expect-error ts-migrate(1238) FIXME: Unable to resolve signature of class decorator whe... Remove this comment to see the full error message
@connect((state) => {
  return {
    genesetsUI: (state as any).genesetsUI,
    colorAccessor: (state as any).colors.colorAccessor,
  };
})
class GenesetMenus extends React.PureComponent<{}, State> {
  constructor(props: {}) {
    super(props);
    this.state = {};
  }

  activateAddGeneToGenesetMode = () => {
    // @ts-expect-error ts-migrate(2339) FIXME: Property 'dispatch' does not exist on type 'Readon... Remove this comment to see the full error message
    const { dispatch, geneset } = this.props;
    dispatch({
      type: "geneset: activate add new genes mode",
      geneset,
    });
  };

  activateEditGenesetNameMode = () => {
    // @ts-expect-error ts-migrate(2339) FIXME: Property 'dispatch' does not exist on type 'Readon... Remove this comment to see the full error message
    const { dispatch, geneset } = this.props;

    dispatch({
      type: "geneset: activate rename geneset mode",
      data: geneset,
    });
  };

  handleColorByEntireGeneset = () => {
    // @ts-expect-error ts-migrate(2339) FIXME: Property 'dispatch' does not exist on type 'Readon... Remove this comment to see the full error message
    const { dispatch, geneset } = this.props;

    dispatch({
      type: "color by geneset mean expression",
      geneset,
    });
  };

  handleDeleteGeneset = () => {
    // @ts-expect-error ts-migrate(2339) FIXME: Property 'dispatch' does not exist on type 'Readon... Remove this comment to see the full error message
    const { dispatch, geneset } = this.props;
    dispatch(actions.genesetDelete(geneset));
  };

  render() {
    // @ts-expect-error ts-migrate(2339) FIXME: Property 'geneset' does not exist on type 'Readonl... Remove this comment to see the full error message
    const { geneset, genesetsEditable, createText, colorAccessor } = this.props;

    const isColorBy = geneset === colorAccessor;

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
                data-testclass="handleAddNewGeneToGeneset"
                data-testid={`${geneset}:add-new-gene-to-geneset`}
                icon={<Icon icon="plus" iconSize={10} />}
                onClick={this.activateAddGeneToGenesetMode}
                small
                minimal
              />
            </Tooltip2>
            {/* @ts-expect-error ts-migrate(2322) FIXME: Type '{ geneset: any; }' is not assignable to type... Remove this comment to see the full error message */}
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
                    data-testclass="handleDeleteGeneset"
                    data-testid={`${geneset}:delete-geneset`}
                    onClick={this.handleDeleteGeneset}
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
              content={`Color by gene set ${geneset} mean`}
              position={Position.BOTTOM}
              hoverOpenDelay={globals.tooltipHoverOpenDelay}
            >
              <AnchorButton
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
