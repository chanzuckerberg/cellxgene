import React from "react";
import { Button, MenuItem } from "@blueprintjs/core";
import { Select } from "@blueprintjs/select";

class DuplicateCategorySelect extends React.PureComponent {
  constructor(props) {
    super(props);
    this.state = {};
  }

  render() {
    const {
      allCategoryNames,
      categoryToDuplicate,
      handleModalDuplicateCategorySelection,
    } = this.props;
    return (
      <div>
        <p>
          Optionally duplicate all labels & cell assignments from existing
          category into new category:
        </p>
        <Select
          items={
            allCategoryNames ||
            [] /* this is a placeholder, could be  a subcomponent to avoid this */
          }
          filterable={false}
          itemRenderer={(d, { handleClick }) => {
            return (
              <MenuItem
                data-testclass="duplicate-category-dropdown-option"
                onClick={handleClick}
                key={d}
                text={d}
              />
            );
          }}
          noResults={<MenuItem disabled text="No results." />}
          onItemSelect={(d) => {
            handleModalDuplicateCategorySelection(d);
          }}
        >
          {/* children become the popover target; render value here */}
          <Button
            data-testid="duplicate-category-dropdown"
            text={categoryToDuplicate || "None (all cells 'unassigned')"}
            rightIcon="double-caret-vertical"
          />
        </Select>
      </div>
    );
  }
}

export default DuplicateCategorySelect;
