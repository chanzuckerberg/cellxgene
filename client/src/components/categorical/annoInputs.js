import React from "react";
import { connect } from "react-redux";
import {
  Button,
  Tooltip,
  InputGroup,
  Dialog,
  MenuItem,
  Classes,
  Colors
} from "@blueprintjs/core";
import { Select } from "@blueprintjs/select";
import fuzzysort from "fuzzysort";

import * as globals from "../../globals";
import sortedCategoryValues from "./util";
import { AnnotationsHelpers } from "../../util/stateManager";

const filterOntology = (query, genes) =>
  /* fires on load, once, and then for each character typed into the input */
  fuzzysort.go(query, genes, {
    limit: 5,
    threshold: -10000 // don't return bad results
  });

@connect(state => ({
  colorAccessor: state.colors.colorAccessor,
  categoricalSelection: state.categoricalSelection,
  annotations: state.annotations,
  universe: state.universe,
  ontology: state.ontology,
  ontologyLoading: state.ontology?.loading
}))
class AnnoInputs extends React.Component {
  constructor(props) {
    super(props);
    this.state = {};
  }

  render() {
    const { useSuggest, handleUserTyping, newCategoryText } = this.props;
    return (
      <div>
        {useSuggest ? (
          <Suggest
            resetOnSelect
            closeOnSelect
            resetOnClose
            itemDisabled={ontologyLoading ? () => true : () => false}
            noResults={<MenuItem disabled text="No matching genes." />}
            onItemSelect={g => {
              /* this happens on 'enter' */
              this.handleClick(g);
            }}
            initialContent={<MenuItem disabled text="Enter a gene…" />}
            inputProps={{ "data-testid": "gene-search" }}
            inputValueRenderer={() => {
              return "";
            }}
            itemListPredicate={filterGenes}
            onActiveItemChange={item => this.setState({ activeItem: item })}
            itemRenderer={renderGene.bind(this)}
            items={
              world && world.varAnnotations
                ? world.varAnnotations.col(varIndexName).asArray()
                : ["No genes"]
            }
            popoverProps={{ minimal: true }}
          />
        ) : (
          <InputGroup
            autoFocus
            value={newCategoryText}
            intent={
              this.categoryNameError(newCategoryText) ? "warning" : "none"
            }
            onChange={e => handleUserTyping(e)}
            leftIcon="tag"
          />
        )}
      </div>
    );
  }
}

export default AnnoInputs;
