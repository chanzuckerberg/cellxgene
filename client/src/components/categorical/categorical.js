// jshint esversion: 6
import React from "react";
import _ from "lodash";
import { connect } from "react-redux";
import * as globals from "../../globals";
import Category from "./category";

/* Cap the max number of displayed categories */
const truncateCategories = options => {
  const numOptions = _.size(options);
  if (numOptions <= globals.maxCategoricalOptionsToDisplay) {
    return options;
  }
  return _(options)
    .map((v, k) => ({ name: k, val: v }))
    .sortBy("val")
    .slice(numOptions - globals.maxCategoricalOptionsToDisplay)
    .transform((r, v) => {
      r[v.name] = v.val;
    }, {})
    .value();
};

@connect(state => ({
  ranges: _.get(state.controls.world, "summary.obs", null),
  categorySelectionLimit: _.get(
    state.config,
    "parameters.max-category-items",
    globals.configDefaults.parameters["max-category-items"]
  )
}))
class Categories extends React.Component {
  render() {
    const { ranges, categorySelectionLimit } = this.props;
    if (!ranges) return null;

    return (
      <div
        style={{
          padding: globals.leftSidebarSectionPadding
        }}
      >
        <p
          style={Object.assign({}, globals.leftSidebarSectionHeading, {
            marginTop: 4
          })}
        >
          Categorical Metadata
        </p>
        {_.map(ranges, (value, key) => {
          const isColorField = key.includes("color") || key.includes("Color");
          const isSelectableCategory =
            value.options &&
            !isColorField &&
            key !== "name" &&
            value.numOptions < categorySelectionLimit;

          if (isSelectableCategory) {
            const categoryOptions = truncateCategories(value.options);
            return (
              <Category
                key={key}
                metadataField={key}
                values={categoryOptions}
                isTruncated={categoryOptions !== value.options}
              />
            );
          }
          return undefined;
        })}
      </div>
    );
  }
}

export default Categories;
