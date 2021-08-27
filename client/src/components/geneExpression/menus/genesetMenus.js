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

import { createColorQuery } from "../../../util/stateManager/colorHelpers";

import * as globals from "../../../globals";
import actions from "../../../actions";
import AddGeneToGenesetDialogue from "./addGeneToGenesetDialogue";

@connect((state) => ({
  annoMatrix: state.annoMatrix,
  schema: state.annoMatrix?.schema,
  genesetsUI: state.genesetsUI,
  colorAccessor: state.colors.colorAccessor,
  genesets: state.genesets.genesets,
}))
class GenesetMenus extends React.PureComponent {
  constructor(props) {
    super(props);
    this.state = {
      dataIsFetching: false,
    };
  }

  fetchOurData = () => {
    /* 
    fetch our data. serves two purposes:  a) preload upon initial render, 
    and b) set busy/loading indicator on any controls that require the data
    to be present, such as the color-by component.
    */
    const { geneset, schema, annoMatrix, genesets } = this.props;
    this.setState({ dataIsFetching: true });
    annoMatrix
      .fetch(
        createColorQuery(
          "color by geneset mean expression",
          geneset,
          schema,
          genesets
        )
      )
      .then(() => this.setState({ dataIsFetching: false }));
  };

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
    this.fetchOurData(); // just in case data was flushed from cache
    dispatch({
      type: "color by geneset mean expression",
      geneset,
    });
  };

  handleDeleteGeneset = () => {
    const { dispatch, geneset } = this.props;
    dispatch(actions.genesetDelete(geneset));
  };

  render() {
    const { geneset, genesetsEditable, createText, colorAccessor } = this.props;
    const { dataIsFetching } = this.state;

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
                loading={dataIsFetching}
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
