import React from "react";
import { connect } from "react-redux";
import { Button, MenuItem } from "@blueprintjs/core";
import { Select } from "@blueprintjs/select";

@connect(state => ({
  colorAccessor: state.colors.colorAccessor,
  categoricalSelection: state.categoricalSelection,
  annotations: state.annotations,
  universe: state.universe,
  ontology: state.ontology,
  ontologyLoading: state.ontology?.loading
}))
class DuplicateCategorySelect extends React.Component {
  constructor(props) {
    super(props);
    this.state = {};
  }

  render() {
    const {
      allCategoryNames,
      categoryToDuplicate,
      handleModalDuplicateCategorySelection
    } = this.props;
    return (
      <Select
        items={
          allCategoryNames ||
          [] /* this is a placeholder, could be  a subcomponent to avoid this */
        }
        filterable={false}
        itemRenderer={(d, { handleClick }) => {
          return <MenuItem onClick={handleClick} key={d} text={d} />;
        }}
        noResults={<MenuItem disabled text="No results." />}
        onItemSelect={d => {
          handleModalDuplicateCategorySelection(d);
        }}
      >
        {/* children become the popover target; render value here */}
        <Button
          text={categoryToDuplicate || "None (all cells 'unassigned')"}
          rightIcon="double-caret-vertical"
        />
      </Select>
    );
  }
}

export default DuplicateCategorySelect;
