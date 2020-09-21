// jshint esversion: 6
/* rc slider https://www.npmjs.com/package/rc-slider */

import React from "react";
import _ from "lodash";
import { connect } from "react-redux";
import { AnchorButton, Icon, Tooltip, Position } from "@blueprintjs/core";
import { FaChevronRight, FaChevronDown } from "react-icons/fa";
import actions from "../../actions";
import Gene from "./gene";
import { memoize } from "../../util/dataframe/util";
import Truncate from "../util/truncate";
// import TestMiniHisto from "./test_miniHisto";
import * as globals from "../../globals";
// import GenesetMenus from "./menus/genesetMenus";

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
      haveFetched: false,
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
        /* Unlike bulk add... with gene sets this test needs to be within the set... */
        // if (upperUserDefinedGenes.get(upperGene) !== undefined) {
        //   return keepAroundErrorToast(`${upperGene} already exists`);
        // }

        const indexOfGene = upperWorldGenes.get(upperGene);

        if (indexOfGene === undefined) {
          console.log("found a gene that doesn't appear to be a valid name");
          // return;
          // keepAroundErrorToast(
          //   `${
          //     genes[upperGenes.get(upperGene)]
          //   } doesn't appear to be a valid gene name.`
          // );
        }
        return dispatch(
          actions.requestUserDefinedGene(worldGenes[indexOfGene])
        );
      })
    ).then(
      () => dispatch({ type: "bulk user defined gene complete" }),
      () => dispatch({ type: "bulk user defined gene error" })
    );

    this.setState({ haveFetched: true });
    return undefined;
  };

  onGenesetMenuClick = () => {
    const { isOpen, haveFetched } = this.state;
    this.setState({ isOpen: !isOpen });
    if (!haveFetched) {
      // this.fetchGenes(); TODO
    }
  };

  onColorChangeClick = () => {
    const { dispatch, setName } = this.props;
    dispatch({
      type: "color by gene set",
      colorAccessor: setName,
    });
  };

  render() {
    const { setName, setGenes, isColorAccessor } = this.props;
    const { isOpen } = this.state;
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
            onKeyPress={/* todo_genesets */ () => {}}
            style={{
              cursor: "pointer",
            }}
            onClick={this.onGenesetMenuClick}
          >
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
            <Truncate>
              <span
                style={{
                  maxWidth:
                    globals.leftSidebarWidth -
                    150 /* todo_genesets this magic number determines how much of a long geneset name we see, and will be tweaked as we build */,
                }}
                data-testid={`${setName}:geneset-label`}
              >
                {setName}
              </span>
            </Truncate>
          </span>
          <div>
            {/* <TestMiniHisto /> */}
            {/* <GenesetMenus genesetsEditable geneset={setName} /> */}
            <Tooltip
              content="Color by geneset"
              position={Position.LEFT}
              usePortal
              hoverOpenDelay={globals.tooltipHoverOpenDelay}
              modifiers={{
                preventOverflow: { enabled: false },
                hide: { enabled: false },
              }}
            >
              <AnchorButton
                disabled
                data-testclass="colorby"
                data-testid={`colorby-${setName}`}
                onClick={this.onColorChangeClick}
                active={isColorAccessor}
                intent={isColorAccessor ? "primary" : "none"}
                icon={<Icon icon="tint" iconSize={16} />}
              />
            </Tooltip>
          </div>
        </div>

        {isOpen
          ? _.map(setGenes, (gene) => {
              return <Gene key={gene} gene={gene} geneset={setName} />;
            })
          : null}
      </div>
    );
  }
}

export default GeneSet;
