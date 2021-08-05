import React from "react";
import { connect } from "react-redux";

import Column from "./column";

@connect((state) => ({
  layoutChoice: state.layoutChoice,
  genesets: state.genesets.genesets,
  dotplot: state.dotplot,
}))
class Dotplot extends React.Component {
  constructor(props) {
    super(props);
    const viewport = this.getViewportDimensions();
    this.dotplotTopPadding = 120;
    this.dotplotLeftPadding = 170;
    this.rowColumnSize = 15;
    this.dotplotBrowserScalingFactor = 0.65;

    this.state = {
      viewport,
    };
  }

  componentDidMount() {
    window.addEventListener("resize", this.handleResize);
    /* chrome only, which is support matrix as of 2021 */
    document.body.style.zoom = `${this.dotplotBrowserScalingFactor * 100}%`;
  }

  componentWillUnmount() {
    window.removeEventListener("resize", this.handleResize);
    document.body.style.zoom = "100%";
  }

  getViewportDimensions = () => {
    const { viewportRef } = this.props;
    return {
      height: viewportRef.clientHeight,
      width: viewportRef.clientWidth,
    };
  };

  render() {
    const { viewport } = this.state;
    const { genesets, dotplot } = this.props;

    let _geneset = null;
    let _genes = null;

    if (dotplot.column) {
      _geneset = genesets.get(dotplot.column);
      _genes = Array.from(_geneset.genes.keys());
    }

    return (
      <div
        id="dotplot-wrapper"
        style={{
          position: "relative",
          top: 0,
          left: 0,
          zIndex: -9999,
        }}
      >
        <svg
          id="dotplot"
          style={{
            position: "absolute",
            top: 0,
            left: 0,
            zIndex: 1,
          }}
          width={viewport.width * (1 + this.dotplotBrowserScalingFactor)}
          height={viewport.height * (1 + this.dotplotBrowserScalingFactor)}
        >
          <g
            id="dotplot_help_text"
            transform={`translate(${this.dotplotLeftPadding},${this.dotplotTopPadding})`}
          >
            <text>
              {!dotplot.row && "Select a row"}{" "}
              {!dotplot.column && "Select a column"}
            </text>
          </g>
          {/* ASYNC HERE */}
          {dotplot.row && dotplot.column && (
            <g
              id="dotplot_interface_margin"
              transform={`translate(${this.dotplotLeftPadding},${this.dotplotTopPadding})`}
            >
              {/* Acaa1b, Mal, Foxq1 ... across the top of the dotplot */}
              <g id="dotplot_column_labels" transform="translate(14,-13)">
                {_genes.map((_geneSymbol, _geneIndexInGeneset) => (
                    <text
                      key={_geneSymbol}
                      x={0}
                      y={0}
                      transform={`translate(${
                        _geneIndexInGeneset * this.rowColumnSize
                      }) rotate(270)`}
                      style={{ fill: "black", font: "12px Roboto Condensed" }}
                    >
                      {_geneSymbol}
                    </text>
                  ))}
              </g>
              {/* loop over genes in the geneset, */}
              <g id="dotplot_columns">
                {_genes.map((_geneSymbol, _geneIndexInGeneset) => (
                    <Column
                      key={_geneSymbol}
                      _geneSymbol={_geneSymbol}
                      _geneIndex={_geneIndexInGeneset}
                      viewport={viewport}
                      rowColumnSize={this.rowColumnSize}
                      metadataField={dotplot.row}
                    />
                  ))}
              </g>
            </g>
          )}
        </svg>
      </div>
    );
  }
}

export default Dotplot;
