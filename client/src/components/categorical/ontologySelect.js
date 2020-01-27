import React from "react";
import { connect } from "react-redux";
import { Button, MenuItem } from "@blueprintjs/core";
import { Select } from "@blueprintjs/select";
import fuzzysort from "fuzzysort";

const filterOntology = (query, ontology) =>
  /* fires on load, once, and then for each character typed into the input */
  fuzzysort.go(query, ontology, {
    limit: 100,
    threshold: -10000 // don't return bad results
  });

@connect()
class ChooseOntologySelect extends React.Component {
  constructor(props) {
    super(props);
    this.state = {};
  }

  render() {
    const {
      ontology,
      categoryToDuplicate,
      handleChooseOntologyTermFromDropdown
    } = this.props;
    return (
      <div>
        <Select
          items={
            ontology?.terms ||
            [] /* this is a placeholder, could be  a subcomponent to avoid this */
          }
          filterable
          itemListPredicate={filterOntology}
          itemRenderer={(d, { handleClick }) => {
            return (
              <MenuItem onClick={handleClick} key={d.target} text={d.target} />
            );
          }}
          noResults={<MenuItem disabled text="No results." />}
          onItemSelect={d => {
            handleChooseOntologyTermFromDropdown(d);
          }}
        >
          {/* children become the popover target; render value here */}
          <Button
            text={categoryToDuplicate || "Choose an Ontology Term"}
            rightIcon="double-caret-vertical"
          />
        </Select>
      </div>
    );
  }
}

export default ChooseOntologySelect;
