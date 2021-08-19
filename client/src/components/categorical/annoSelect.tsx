import React from "react";
import { Button, MenuItem } from "@blueprintjs/core";
import { Select } from "@blueprintjs/select";

type State = any;

class DuplicateCategorySelect extends React.PureComponent<{}, State> {
  constructor(props: {}) {
    super(props);
    this.state = {};
  }

  render() {
    const {
      // @ts-expect-error ts-migrate(2339) FIXME: Property 'allCategoryNames' does not exist on type... Remove this comment to see the full error message
      allCategoryNames,
      // @ts-expect-error ts-migrate(2339) FIXME: Property 'categoryToDuplicate' does not exist on t... Remove this comment to see the full error message
      categoryToDuplicate,
      // @ts-expect-error ts-migrate(2339) FIXME: Property 'handleModalDuplicateCategorySelection' d... Remove this comment to see the full error message
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
                // @ts-expect-error ts-migrate(2322) FIXME: Type 'unknown' is not assignable to type 'Key | nu... Remove this comment to see the full error message
                key={d}
                // @ts-expect-error ts-migrate(2322) FIXME: Type 'unknown' is not assignable to type 'ReactNod... Remove this comment to see the full error message
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
