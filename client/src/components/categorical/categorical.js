// jshint esversion: 6
import React from "react";
import _ from "lodash";
import { connect } from "react-redux";
import * as globals from "../../globals";

import Category from "./category";

@connect(state => {
  const ranges = _.get(state.controls.world, "summary.obs", null);
  return {
    ranges
  };
})
class Categories extends React.Component {
  constructor(props) {
    super(props);
    this.state = {};
  }

  render() {
    const { ranges } = this.props;
    if (!ranges) return null;

    return (
      <div style={{}}>
        <p
          style={Object.assign({}, globals.leftSidebarSectionHeading, {
            marginTop: 14
          })}
        >
          Categorical Metadata
        </p>
        {_.map(ranges, (value, key) => {
          const isColorField = key.includes("color") || key.includes("Color");
          if (value.options && !isColorField && key !== "name") {
            return (
              <Category key={key} metadataField={key} values={value.options} />
            );
          }
          return undefined;
        })}
      </div>
    );
  }
}

export default Categories;
