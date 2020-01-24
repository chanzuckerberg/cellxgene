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
    handleCreateArbitraryLabel,
    handleTextChange,
    isTextInvalid,
    isTextInvalidErrorMessage
  } = props;
  return (
    <Suggest
      fill
      resetOnSelect
      closeOnSelect
      resetOnClose
      createNewItemFromQuery={() => {}}
      createNewItemRenderer={userInputStr => {
        return isTextInvalid?.(userInputStr) ? (
          <MenuItem disabled text={isTextInvalidErrorMessage(userInputStr)} />
        ) : (
          <MenuItem
            icon="add"
            text="Create a label not in the ontology"
            active
            onClick={() => {
              handleCreateArbitraryLabel(userInputStr);
            }}
            shouldDismissPopover={false}
          />
        );
      }}
      itemDisabled={ontologyLoading ? () => true : () => false}
      noResults={<MenuItem disabled text="No matching ontology identifier" />}
      onItemSelect={handleChoice}
      initialContent={
        <MenuItem disabled text="Enter an ontology identifier…" />
      }
      inputProps={{ "data-testid": "gene-search" }}
      inputValueRenderer={t => t.target}
      itemListPredicate={filterOntology}
      onActiveItemChange={handleItemChange}
      itemRenderer={renderListItem.bind(this)}
      items={!ontologyLoading && ontology ? ontology : ["No ontology loaded"]}
      popoverProps={{ minimal: true }}
      onQueryChange={(s, e) => {
        // undefined event means resetOnSelect
        if (e !== undefined) handleTextChange?.(s);
      }}
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
      onChange={e => handleTextChange?.(e.target.value)}
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
    const { handleTextChange, text } = this.props;
    return (
      <div>
        <VanillaInput text={text} handleTextChange={handleTextChange} />
      </div>
    );
  }
}

export default AnnoInputs;
