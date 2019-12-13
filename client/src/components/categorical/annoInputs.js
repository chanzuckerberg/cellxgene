import React from "react";
import { connect } from "react-redux";
import { InputGroup, MenuItem } from "@blueprintjs/core";
import { Suggest } from "@blueprintjs/select";
import fuzzysort from "fuzzysort";

const filterOntology = (query, genes) =>
  /* fires on load, once, and then for each character typed into the input */
  fuzzysort.go(query, genes, {
    limit: 5,
    threshold: -10000 // don't return bad results
  });

const renderListItem = (fuzzySortResult, { handleClick, modifiers }) => {
  if (!modifiers.matchesPredicate) {
    return null;
  }
  /* the fuzzysort wraps the object with other properties, like a score */
  const geneName = fuzzySortResult.target;

  return (
    <MenuItem
      active={modifiers.active}
      disabled={modifiers.disabled}
      data-testid={`suggest-menu-item-${geneName}`}
      // Use of annotations in this way is incorrect and dataset specific.
      // See https://github.com/chanzuckerberg/cellxgene/issues/483
      // label={gene.n_counts}
      key={geneName}
      onClick={g =>
        /* this fires when user clicks a menu item */
        handleClick(g)
      }
      text={geneName}
    />
  );
};

const AnnoSuggest = (world, varIndexName, ontologyLoading) => {
  return (
    <Suggest
      resetOnSelect
      closeOnSelect
      resetOnClose
      itemDisabled={ontologyLoading ? () => true : () => false}
      noResults={<MenuItem disabled text="No matching genes." />}
      onItemSelect={g => {
        /* this happens on 'enter' */
        console.log("annoInput Suggest sees an item selected...", g);
      }}
      initialContent={<MenuItem disabled text="Enter a gene…" />}
      inputProps={{ "data-testid": "gene-search" }}
      inputValueRenderer={() => {
        return "";
      }}
      itemListPredicate={filterOntology}
      onActiveItemChange={item => this.setState({ activeItem: item })}
      itemRenderer={renderListItem.bind(this)}
      items={
        world && world.varAnnotations
          ? world.varAnnotations.col(varIndexName).asArray()
          : ["No genes"]
      }
      popoverProps={{ minimal: true }}
    />
  );
};

const VanillaInput = () => {
  return (
    <InputGroup
      autoFocus
      value={newCategoryText}
      intent={this.categoryNameError(newCategoryText) ? "warning" : "none"}
      onChange={e => handleUserTyping(e)}
      leftIcon="tag"
    />
  );
};

@connect(state => ({
  colorAccessor: state.colors.colorAccessor,
  categoricalSelection: state.categoricalSelection,
  annotations: state.annotations,
  universe: state.universe,
  world: state.world,
  ontology: state.ontology,
  ontologyLoading: state.ontology?.loading
}))
class AnnoInputs extends React.Component {
  constructor(props) {
    super(props);
    this.state = {};
  }

  render() {
    const {
      useSuggest,
      handleUserTyping,
      newCategoryText,
      ontologyLoading,
      world
    } = this.props;
    const varIndexName = world.schema.annotations.var.index;

    return (
      <div>
        {useSuggest ? (
          <AnnoSuggest
            world={world}
            varIndexName={varIndexName}
            ontologyLoading={ontologyLoading}
          />
        ) : (
          <VanillaInput />
        )}
      </div>
    );
  }
}

export default AnnoInputs;
