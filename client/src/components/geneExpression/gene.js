// jshint esversion: 6
/* rc slider https://www.npmjs.com/package/rc-slider */

import React from "react";
import { connect } from "react-redux";
// import { FaChevronRight, FaChevronDown } from "react-icons/fa";

import { AnchorButton, Icon } from "@blueprintjs/core";
import Truncate from "../util/truncate";
import HistogramBrush from "../brushableHistogram";

// import TestMiniHisto from "./test_miniHisto";
import * as globals from "../../globals";
import GeneMenus from "./menus/geneMenus";

@connect(() => {
  return {};
})
class Gene extends React.Component {
  constructor(props) {
    super(props);
    this.state = {
      geneIsExpanded: false,
    };
  }

  onColorChangeClick = () => {
    const { dispatch, gene } = this.props;
    dispatch({
      type: "color by expression",
      colorAccessor: gene,
    });
  };

  handleGeneExpandClick = () => {
    const { geneIsExpanded } = this.state;
    this.setState({ geneIsExpanded: !geneIsExpanded });
  };

  render() {
    const { gene } = this.props;
    const { geneIsExpanded } = this.state;
    return (
      <div
        style={{
          display: "flex",
          justifyContent: "space-between",
          alignItems: "baseline",
          marginLeft: 15,
        }}
      >
        <div>
          <Icon
            icon="drag-handle-horizontal"
            iconSize={12}
            style={{
              marginRight: 10,
              cursor: "grab",
              position: "relative",
              top: -2,
            }}
          />
          <span
            role="menuitem"
            tabIndex="0"
            data-testclass="gene-expand"
            data-testid={`${gene}:gene-expand`}
            onKeyPress={/* todo_genesets */ () => {}}
            style={{
              cursor: "pointer",
            }}
            onClick={this.handleGeneExpandClick}
          >
            <Truncate>
              <span
                style={{
                  width:
                    globals.leftSidebarWidth -
                    310 /* todo_genesets this magic number determines how much of a long geneset name we see, and will be tweaked as we build */,
                }}
                data-testid={`${gene}:gene-label`}
              >
                {gene}
              </span>
            </Truncate>
            <HistogramBrush isUserDefined field={gene} mini={!geneIsExpanded} />

            {/* <TestMiniHisto /> */}
          </span>
        </div>

        <div>
          <GeneMenus genesetsEditable gene={gene} />
          <AnchorButton
            minimal
            small
            data-testclass="colorby"
            data-testid={`colorby-${gene}`}
            onClick={/* todo gene sets */ () => {}}
            active={false /* todo gene sets */}
            intent="none"
            icon={<Icon icon="tint" iconSize={12} />}
          />
          <AnchorButton
            minimal
            small
            data-testclass="maximize"
            data-testid={`maximize-${gene}`}
            onClick={/* todo gene sets */ () => {}}
            active={false /* todo gene sets */}
            intent="none"
            icon={<Icon icon="maximize" iconSize={10} />}
          />
        </div>
      </div>
    );
  }
}

export default Gene;

// {geneIsExpanded ? (
//   <FaChevronDown
//     data-testclass="gene-expand-is-expanded"
//     style={{ fontSize: 10, marginRight: 5 }}
//   />
// ) : (
//   <FaChevronRight
//     data-testclass="gene-expand-is-not-expanded"
//     style={{ fontSize: 10, marginRight: 5 }}
//   />
// )}
