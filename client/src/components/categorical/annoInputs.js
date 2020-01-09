import React from "react";
import { connect } from "react-redux";
import { InputGroup, MenuItem } from "@blueprintjs/core";
import { Suggest } from "@blueprintjs/select";
import fuzzysort from "fuzzysort";

const filterOntology = (query, ontology) =>
  /* fires on load, once, and then for each character typed into the input */
  fuzzysort.go(query, ontology, {
    limit: 100,
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

const AnnoSuggest = props => {
  const {
    ontologyLoading,
    handleItemChange,
    handleChoice,
    ontology,
    text
  } = props;
  return (
    <Suggest
      fill
      autoFocus
      resetOnSelect
      closeOnSelect
      resetOnClose
      itemDisabled={ontologyLoading ? () => true : () => false}
      noResults={<MenuItem disabled text="No matching ontology identifier" />}
      onItemSelect={handleChoice}
      initialContent={
        <MenuItem disabled text="Enter an ontology identifier…" />
      }
      inputProps={{ "data-testid": "gene-search" }}
      inputValueRenderer={() => {
        return text || "";
      }}
      itemListPredicate={filterOntology}
      onActiveItemChange={handleItemChange}
      itemRenderer={renderListItem.bind(this)}
      items={!ontologyLoading && ontology ? ontology : ["No ontology loaded"]}
      popoverProps={{ minimal: true }}
    />
  );
};

const VanillaInput = props => {
  const { text, handleTextChange } = props;
  return (
    <InputGroup
      autoFocus
      value={text}
      intent="none"
      onChange={e => handleTextChange(e)}
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
  ontology: state.ontology?.terms,
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
      handleTextChange,
      text,
      ontologyLoading,
      world,
      handleChoice,
      handleItemChange,
      ontology
    } = this.props;
    return (
      <div>
        {/* Le sigh https://github.com/palantir/blueprint/issues/2864 */}
        {useSuggest ? (
          <AnnoSuggest
            world={world}
            text={text}
            handleChoice={handleChoice}
            ontologyLoading={ontologyLoading}
            handleItemChange={handleItemChange}
            ontology={ontology}
          />
        ) : (
          <VanillaInput text={text} handleTextChange={handleTextChange} />
        )}
      </div>
    );
  }
}

export default AnnoInputs;
