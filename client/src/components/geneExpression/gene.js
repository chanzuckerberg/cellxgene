// jshint esversion: 6
/* rc slider https://www.npmjs.com/package/rc-slider */

import React from "react";
import { connect } from "react-redux";
import { FaChevronRight, FaChevronDown } from "react-icons/fa";

import { AnchorButton } from "@blueprintjs/core";
import Truncate from "../util/truncate";

import TestMiniHisto from "./test_miniHisto";
import * as globals from "../../globals";

@connect(() => {
  return {};
})
class Gene extends React.Component {
  constructor(props) {
    super(props);
    this.state = {
      geneIsExpanded: false,
      haveFetched: false,
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
    const { geneIsExpanded, haveFetched } = this.state;
    this.setState({ geneIsExpanded: !geneIsExpanded });
    if (!haveFetched) {
      // this.fetchGene(); TODO
    }
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
          {geneIsExpanded ? (
            <FaChevronDown
              data-testclass="gene-expand-is-expanded"
              style={{ fontSize: 10, marginRight: 5 }}
            />
          ) : (
            <FaChevronRight
              data-testclass="gene-expand-is-not-expanded"
              style={{ fontSize: 10, marginRight: 5 }}
            />
          )}
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
          <TestMiniHisto />
        </span>

        <AnchorButton
          minimal
          data-testclass="colorby"
          data-testid={`colorby-${gene}`}
          onClick={/* todo gene sets */ () => {}}
          active={false /* todo gene sets */}
          intent="none"
          icon="tint"
        />
      </div>
    );
  }
}

export default Gene;
