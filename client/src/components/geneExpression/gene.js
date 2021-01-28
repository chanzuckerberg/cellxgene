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
import actions from "../../actions";

// import GeneMenus from "./menus/geneMenus";

@connect((state, ownProps) => {
  const { gene } = ownProps;

  return {
    isColorAccessor: state.colors.colorAccessor === gene,
    isScatterplotXXaccessor: state.controls.scatterplotXXaccessor === gene,
    isScatterplotYYaccessor: state.controls.scatterplotYYaccessor === gene,
  };
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
    dispatch(actions.requestSingleGeneExpressionCountsForColoringPOST(gene));
  };

  handleGeneExpandClick = () => {
    const { geneIsExpanded } = this.state;
    this.setState({ geneIsExpanded: !geneIsExpanded });
  };

  handleSetGeneAsScatterplotX = () => {
    const { dispatch, gene } = this.props;
    dispatch({
      type: "set scatterplot x",
      data: gene,
    });
  };

  handleSetGeneAsScatterplotY = () => {
    const { dispatch, gene } = this.props;
    dispatch({
      type: "set scatterplot y",
      data: gene,
    });
  };

  handleDeleteGeneFromSet = () => {
    const { dispatch, gene, geneset } = this.props;

    dispatch({
      type: "geneset: del genes",
      name: geneset,
      genes: [gene],
    });
  };

  render() {
    const {
      gene,
      isColorAccessor,
      isScatterplotXXaccessor,
      isScatterplotYYaccessor,
      isDiffexp,
    } = this.props;
    const { geneIsExpanded } = this.state;
    return (
      <div>
        <div
          style={{
            marginLeft: 5,
            marginRight: 0,
            marginTop: 2,
            display: "flex",
            justifyContent: "space-between",
            alignItems: "center",
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
              display: "flex",
              justifyContent: "space-between",
            }}
          >
            <Icon
              icon="drag-handle-horizontal"
              iconSize={12}
              style={{
                marginRight: 7,
                cursor: "grab",
                position: "relative",
                top: 3,
              }}
            />
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
            {!geneIsExpanded ? (
              <HistogramBrush isUserDefined field={gene} mini />
            ) : null}
          </span>
          <span>
            {!isDiffexp ? (
              <AnchorButton
                minimal
                small
                data-testid={`delete-from-geneset-${gene}`}
                onClick={this.handleDeleteGeneFromSet}
                intent="danger"
                style={{ fontWeight: 700, marginRight: 2 }}
                icon={<Icon icon="trash" iconSize={10} />}
              />
            ) : null}
            <AnchorButton
              minimal
              small
              data-testid={`plot-x-${gene}`}
              onClick={this.handleSetGeneAsScatterplotX}
              active={isScatterplotXXaccessor}
              intent={isScatterplotXXaccessor ? "primary" : "none"}
              style={{ fontWeight: 700, marginRight: 2 }}
            >
              x
            </AnchorButton>
            <AnchorButton
              minimal
              small
              data-testid={`plot-y-${gene}`}
              onClick={this.handleSetGeneAsScatterplotY}
              active={isScatterplotYYaccessor}
              intent={isScatterplotYYaccessor ? "primary" : "none"}
              style={{ fontWeight: 700, marginRight: 2 }}
            >
              y
            </AnchorButton>
            <AnchorButton
              minimal
              small
              data-testclass="maximize"
              data-testid={`maximize-${gene}`}
              onClick={this.handleGeneExpandClick}
              active={geneIsExpanded}
              intent="none"
              icon={<Icon icon="maximize" iconSize={10} />}
              style={{ marginRight: 2 }}
            />
            <AnchorButton
              minimal
              small
              data-testclass="colorby"
              data-testid={`colorby-${gene}`}
              onClick={this.onColorChangeClick}
              active={isColorAccessor}
              intent={isColorAccessor ? "primary" : "none"}
              icon={<Icon icon="tint" iconSize={12} />}
            />
          </span>
          {/* <TestMiniHisto /> */}
        </div>
        {geneIsExpanded ? <HistogramBrush isUserDefined field={gene} /> : null}
        {/* <GeneMenus genesetsEditable gene={gene} /> */}
      </div>
    );
  }
}

export default Gene;
