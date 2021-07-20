import React from "react";
import { FaChevronRight, FaChevronDown } from "react-icons/fa";
import Gene from "./gene";
import Truncate from "../util/truncate";
import * as globals from "../../globals";
import GenesetMenus from "./menus/genesetMenus";
import EditGenesetNameDialogue from "./menus/editGenesetNameDialogue";
import HistogramBrush from "../brushableHistogram";

import { diffexpPopNamePrefix1, diffexpPopNamePrefix2 } from "../../globals";

class GeneSet extends React.Component {
  constructor(props) {
    super(props);
    this.state = {
      isOpen: false,
    };
  }

  onGenesetMenuClick = () => {
    const { isOpen } = this.state;
    this.setState({ isOpen: !isOpen });
  };

  renderGenes() {
    const { setName, setGenes } = this.props;
    const setGenesNames = [...setGenes.keys()];
    return (
      <div data-testclass="gene-set-genes">
        {setGenesNames.map((gene) => {
          const { geneDescription } = setGenes.get(gene);
          return (
            <Gene
              key={gene}
              gene={gene}
              geneDescription={geneDescription}
              geneset={setName}
            />
          );
        })}
      </div>
    );
  }

  render() {
    const { setName, genesetDescription, setGenes } = this.props;
    const { isOpen } = this.state;
    const genesetNameLengthVisible = 150; /* this magic number determines how much of a long geneset name we see */
    const genesetIsEmpty = setGenes.size === 0;
    let testClass = "geneset-expand";

    if (setName.includes(diffexpPopNamePrefix1))
      testClass = "pop-1-geneset-expand";
    else if (setName.includes(diffexpPopNamePrefix2))
      testClass = "pop-2-geneset-expand";

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
            data-testclass={testClass}
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
            <Truncate
              tooltipAddendum={
                genesetDescription ? `: ${genesetDescription}` : ""
              }
            >
              <span
                style={{
                  maxWidth: globals.leftSidebarWidth - genesetNameLengthVisible,
                }}
                data-testid={`${setName}:geneset-name`}
              >
                {setName}
              </span>
            </Truncate>
            {isOpen ? (
              <FaChevronDown
                data-testclass="geneset-expand-is-expanded"
                style={{ fontSize: 10, marginLeft: 5 }}
              />
            ) : (
              <FaChevronRight
                data-testclass="geneset-expand-is-not-expanded"
                style={{ fontSize: 10, marginLeft: 5 }}
              />
            )}
          </span>
          <div>
            <GenesetMenus isOpen={isOpen} genesetsEditable geneset={setName} />
          </div>
        </div>

        <div style={{ marginLeft: 15, marginTop: 5, marginRight: 0 }}>
          {isOpen && genesetIsEmpty && (
            <p style={{ fontStyle: "italic", color: "lightgrey" }}>
              No genes to display
            </p>
          )}
        </div>
        {isOpen && !genesetIsEmpty && (
          <HistogramBrush
            isGeneSetSummary
            field={setName}
            setGenes={setGenes}
          />
        )}
        {isOpen && !genesetIsEmpty && this.renderGenes()}
        <EditGenesetNameDialogue
          parentGeneset={setName}
          parentGenesetDescription={genesetDescription}
        />
      </div>
    );
  }
}

export default GeneSet;
