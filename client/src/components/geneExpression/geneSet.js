import React from "react";
import _ from "lodash";
import { connect } from "react-redux";
import { Tooltip, Position, Switch } from "@blueprintjs/core";
import { FaChevronRight, FaChevronDown } from "react-icons/fa";
import actions from "../../actions";
import Gene from "./gene";
import { memoize } from "../../util/dataframe/util";
import Truncate from "../util/truncate";
import * as globals from "../../globals";
import GenesetMenus from "./menus/genesetMenus";

@connect((state, ownProps) => {
  return {
    world: state.world,
    userDefinedGenes: state.controls.userDefinedGenes,
    userDefinedGenesLoading: state.controls.userDefinedGenesLoading,
    isColorAccessor: state.colors.colorAccessor === ownProps.setName,
  };
})
class GeneSet extends React.Component {
  _memoGenesToUpper = memoize(this._genesToUpper, (arr) => arr);

  constructor(props) {
    super(props);
    this.state = {
      isOpen: false,
      toggleSummaryHisto: false,
    };
  }

  _genesToUpper = (listGenes) => {
    // Has to be a Map to preserve index
    const upperGenes = new Map();
    for (let i = 0, { length } = listGenes; i < length; i += 1) {
      upperGenes.set(listGenes[i].toUpperCase(), i);
    }

    return upperGenes;
  };

  fetchGenes = () => {
    const { world, dispatch, setGenes } = this.props;
    const varIndexName = world.schema.annotations.var.index;

    const worldGenes = world.varAnnotations.col(varIndexName).asArray();

    const upperGenes = this._genesToUpper(setGenes);
    const upperWorldGenes = this._memoGenesToUpper(worldGenes);

    dispatch({ type: "bulk user defined gene start" });

    Promise.all(
      [...upperGenes.keys()].map((upperGene) => {
        const indexOfGene = upperWorldGenes.get(upperGene);

        if (indexOfGene === undefined) {
          console.log("found a gene that doesn't appear to be a valid name");
        }
        return dispatch(
          actions.requestUserDefinedGene(worldGenes[indexOfGene])
        );
      })
    ).then(
      () => dispatch({ type: "bulk user defined gene complete" }),
      () => dispatch({ type: "bulk user defined gene error" })
    );

    return undefined;
  };

  onGenesetMenuClick = () => {
    const { isOpen } = this.state;
    this.setState({ isOpen: !isOpen });
  };

  onColorChangeClick = () => {
    //   const { dispatch, setName } = this.props;
    //   dispatch({
    //     type: "color by gene set",
    //     colorAccessor: setName,
    //   });
  };

  toggleSummaryHisto = () => {
    const { toggleSummaryHisto } = this.state;
    this.setState({ toggleSummaryHisto: !toggleSummaryHisto });
  };

  render() {
    const { setName, setGenes, isDiffexp } = this.props;
    const { isOpen, toggleSummaryHisto } = this.state;
    const genesetNameLengthVisible = 120; /* this magic number determines how much of a long geneset name we see */
    return (
      <div style={{ marginBottom: 3 }}>
        <div
          style={{
            display: "flex",
            justifyContent: "space-between",
            alignItems: "baseline",
          }}
        >
          <span
            role="menuitem"
            tabIndex="0"
            data-testclass="geneset-expand"
            data-testid={`${setName}:geneset-expand`}
            onKeyPress={
              /* TODO(colinmegill): #2101: click handler on span */ () => {}
            }
            style={{
              cursor: "pointer",
              userSelect: "none",
            }}
            onClick={this.onGenesetMenuClick}
          >
            <Truncate>
              <span
                style={{
                  maxWidth: globals.leftSidebarWidth - genesetNameLengthVisible,
                }}
                data-testid={`${setName}:geneset-label`}
              >
                {setName}
              </span>
            </Truncate>
            {isOpen ? (
              <FaChevronDown
                data-testclass="geneset-expand-is-expanded"
                style={{ fontSize: 10, marginRight: 5 }}
              />
            ) : (
              <FaChevronRight
                data-testclass="geneset-expand-is-not-expanded"
                style={{ fontSize: 10, marginRight: 5 }}
              />
            )}
          </span>
          <div>
            <GenesetMenus genesetsEditable geneset={setName} />
          </div>
        </div>

        <div style={{ marginLeft: 15, marginTop: 5, marginRight: 0 }}>
          {isOpen ? (
            <Tooltip
              content="Aggregate all genes in this geneset, and allow coloring by and selecting on the entire set."
              position={Position.BOTTOM_RIGHT}
              usePortal
              hoverOpenDelay={globals.tooltipHoverOpenDelay}
              modifiers={{
                preventOverflow: { enabled: false },
                hide: { enabled: false },
              }}
            >
              <Switch
                style={{ fontStyle: "italic" }}
                checked={toggleSummaryHisto}
                label="Show mean expression across entire set"
                alignIndicator="left"
                innerLabel="off"
                innerLabelChecked="on"
                onChange={this.toggleSummaryHisto}
              />
            </Tooltip>
          ) : null}
        </div>

        {isOpen && !toggleSummaryHisto
          ? _.map(setGenes, (gene) => {
              return (
                <Gene
                  key={gene}
                  gene={gene}
                  geneset={setName}
                  isDiffexp={isDiffexp}
                />
              );
            })
          : null}
      </div>
    );
  }
}

export default GeneSet;
