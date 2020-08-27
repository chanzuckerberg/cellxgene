import React from "react";
import { connect } from "react-redux";
import {
  Button,
  Menu,
  MenuItem,
  Popover,
  Position,
  Icon,
  PopoverInteractionKind,
} from "@blueprintjs/core";

@connect(() => ({}))
class GenesetMenus extends React.PureComponent {
  constructor(props) {
    super(props);
    this.state = {};
  }

  handleDeleteGeneFromGeneset = () => {
    // const { dispatch, geneset } = this.props;
    // dispatch(actions.geneDeleteFromGenesetAction(geneset, gene));
    // todo genesets
  };

  setGeneScatterplotX = () => {
    // todo genesets
  };

  setGeneScatterplotY = () => {
    // todo genesets
  };

  render() {
    const { geneset, genesetsEditable, gene } = this.props;

    return (
      <>
        {genesetsEditable ? (
          <>
            <Popover
              interactionKind={PopoverInteractionKind.HOVER}
              boundary="window"
              position={Position.BOTTOM}
              content={
                <Menu>
                  <MenuItem
                    icon={<Icon icon="scatter-plot" iconSize={16} />}
                    intent="none"
                    data-testclass="setGeneScatterplotX"
                    data-testid={`${geneset}:set-gene-as-scatterplot-x`}
                    onClick={this.setGeneScatterplotX}
                    text={`X: Set ${gene}'s values as as 'x' on a scatterplot`}
                  />
                  <MenuItem
                    icon={<Icon icon="scatter-plot" iconSize={16} />}
                    intent="none"
                    data-testclass="setGeneScatterplotY"
                    data-testid={`${geneset}:set-gene-as-scatterplot-y`}
                    onClick={this.setGeneScatterplotY}
                    text={`Y: Set ${gene}'s values as as 'y' on a scatterplot`}
                  />
                  <MenuItem
                    icon="delete"
                    intent="danger"
                    data-testclass="handleDeleteGeneFromGeneset"
                    data-testid={`${geneset}:delete-gene-from-geneset`}
                    onClick={this.handleDeleteGeneFromGeneset}
                    text="Remove this gene from the geneset"
                  />
                </Menu>
              }
            >
              <Button
                style={{ marginLeft: 0, marginRight: 5 }}
                data-testclass="seeActions"
                data-testid={`${gene}:see-actions`}
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

export default GenesetMenus;
