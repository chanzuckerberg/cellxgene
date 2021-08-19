import React from "react";
import { connect } from "react-redux";

import { Button, Icon } from "@blueprintjs/core";
import Truncate from "../util/truncate";
import HistogramBrush from "../brushableHistogram";

import actions from "../../actions";

const MINI_HISTOGRAM_WIDTH = 110;

type State = any;

// @ts-expect-error ts-migrate(1238) FIXME: Unable to resolve signature of class decorator whe... Remove this comment to see the full error message
@connect((state, ownProps) => {
  // @ts-expect-error ts-migrate(2339) FIXME: Property 'gene' does not exist on type '{}'.
  const { gene } = ownProps;

  return {
    isColorAccessor: (state as any).colors.colorAccessor === gene,
    isScatterplotXXaccessor:
      (state as any).controls.scatterplotXXaccessor === gene,
    isScatterplotYYaccessor:
      (state as any).controls.scatterplotYYaccessor === gene,
  };
})
class Gene extends React.Component<{}, State> {
  constructor(props: {}) {
    super(props);
    this.state = {
      geneIsExpanded: false,
    };
  }

  onColorChangeClick = () => {
    // @ts-expect-error ts-migrate(2339) FIXME: Property 'dispatch' does not exist on type 'Readon... Remove this comment to see the full error message
    const { dispatch, gene } = this.props;
    dispatch(actions.requestSingleGeneExpressionCountsForColoringPOST(gene));
  };

  handleGeneExpandClick = () => {
    const { geneIsExpanded } = this.state;
    this.setState({ geneIsExpanded: !geneIsExpanded });
  };

  handleSetGeneAsScatterplotX = () => {
    // @ts-expect-error ts-migrate(2339) FIXME: Property 'dispatch' does not exist on type 'Readon... Remove this comment to see the full error message
    const { dispatch, gene } = this.props;
    dispatch({
      type: "set scatterplot x",
      data: gene,
    });
  };

  handleSetGeneAsScatterplotY = () => {
    // @ts-expect-error ts-migrate(2339) FIXME: Property 'dispatch' does not exist on type 'Readon... Remove this comment to see the full error message
    const { dispatch, gene } = this.props;
    dispatch({
      type: "set scatterplot y",
      data: gene,
    });
  };

  handleDeleteGeneFromSet = () => {
    // @ts-expect-error ts-migrate(2339) FIXME: Property 'dispatch' does not exist on type 'Readon... Remove this comment to see the full error message
    const { dispatch, gene, geneset } = this.props;
    dispatch(actions.genesetDeleteGenes(geneset, [gene]));
  };

  render() {
    const {
      // @ts-expect-error ts-migrate(2339) FIXME: Property 'gene' does not exist on type 'Readonly<{... Remove this comment to see the full error message
      gene,
      // @ts-expect-error ts-migrate(2339) FIXME: Property 'geneDescription' does not exist on type ... Remove this comment to see the full error message
      geneDescription,
      // @ts-expect-error ts-migrate(2339) FIXME: Property 'isColorAccessor' does not exist on type ... Remove this comment to see the full error message
      isColorAccessor,
      // @ts-expect-error ts-migrate(2339) FIXME: Property 'isScatterplotXXaccessor' does not exist ... Remove this comment to see the full error message
      isScatterplotXXaccessor,
      // @ts-expect-error ts-migrate(2339) FIXME: Property 'isScatterplotYYaccessor' does not exist ... Remove this comment to see the full error message
      isScatterplotYYaccessor,
      // @ts-expect-error ts-migrate(2339) FIXME: Property 'quickGene' does not exist on type 'Reado... Remove this comment to see the full error message
      quickGene,
      // @ts-expect-error ts-migrate(2339) FIXME: Property 'removeGene' does not exist on type 'Read... Remove this comment to see the full error message
      removeGene,
    } = this.props;
    const { geneIsExpanded } = this.state;
    const geneSymbolWidth = 60 + (geneIsExpanded ? MINI_HISTOGRAM_WIDTH : 0);

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
          <div
            role="menuitem"
            // @ts-expect-error ts-migrate(2322) FIXME: Type 'string' is not assignable to type 'number | ... Remove this comment to see the full error message
            tabIndex="0"
            data-testclass="gene-expand"
            data-testid={`${gene}:gene-expand`}
            onKeyPress={() => {}}
            style={{
              cursor: "pointer",
              display: "flex",
              justifyContent: "space-between",
              width: "100%",
            }}
          >
            <div>
              {!quickGene && (
                <Icon
                  icon="drag-handle-horizontal"
                  iconSize={12}
                  style={{
                    marginRight: 7,
                    cursor: "grab",
                    position: "relative",
                    top: -1,
                  }}
                />
              )}
              <Truncate
                tooltipAddendum={geneDescription && `: ${geneDescription}`}
              >
                <span
                  style={{
                    width: geneSymbolWidth,
                    display: "inline-block",
                  }}
                  data-testid={`${gene}:gene-label`}
                >
                  {gene}
                </span>
              </Truncate>
            </div>
            {!geneIsExpanded ? (
              <HistogramBrush
                // @ts-expect-error ts-migrate(2769) FIXME: No overload matches this call.
                isUserDefined
                field={gene}
                mini
                width={MINI_HISTOGRAM_WIDTH}
              />
            ) : null}
          </div>
          <div style={{ flexShrink: 0, marginLeft: 2 }}>
            <Button
              minimal
              small
              data-testid={`delete-from-geneset:${gene}`}
              onClick={
                quickGene ? removeGene(gene) : this.handleDeleteGeneFromSet
              }
              intent="none"
              style={{ fontWeight: 700, marginRight: 2 }}
              icon={<Icon icon="trash" iconSize={10} />}
            />
            <Button
              minimal
              small
              data-testid={`plot-x-${gene}`}
              onClick={this.handleSetGeneAsScatterplotX}
              active={isScatterplotXXaccessor}
              intent={isScatterplotXXaccessor ? "primary" : "none"}
              style={{ fontWeight: 700, marginRight: 2 }}
            >
              x
            </Button>
            <Button
              minimal
              small
              data-testid={`plot-y-${gene}`}
              onClick={this.handleSetGeneAsScatterplotY}
              active={isScatterplotYYaccessor}
              intent={isScatterplotYYaccessor ? "primary" : "none"}
              style={{ fontWeight: 700, marginRight: 2 }}
            >
              y
            </Button>
            <Button
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
            <Button
              minimal
              small
              data-testclass="colorby"
              data-testid={`colorby-${gene}`}
              onClick={this.onColorChangeClick}
              active={isColorAccessor}
              intent={isColorAccessor ? "primary" : "none"}
              icon={<Icon icon="tint" iconSize={12} />}
            />
          </div>
        </div>
        {/* @ts-expect-error ts-migrate(2769) FIXME: No overload matches this call. */}
        {geneIsExpanded && <HistogramBrush isUserDefined field={gene} />}
      </div>
    );
  }
}

export default Gene;
