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
  ranges: _.get(state.controls.world, "summary.obs", null)
}))
class Categories extends React.Component {
  render() {
    const { ranges } = this.props;
    if (!ranges) return null;

    return (
      <div
        style={{
          width: 310,
          marginRight: 40,
          paddingRight: 20,
          flexShrink: 0
        }}
      >
        <p> Categorical Metadata </p>
        {_.map(ranges, (value, key) => {
          const isColorField = key.includes("color") || key.includes("Color");
          if (value.options && !isColorField && key !== "name") {
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
