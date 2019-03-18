// jshint esversion: 6
import React from "react";
import _ from "lodash";
import { connect } from "react-redux";
import * as globals from "../../globals";
import Category from "./category";

@connect(state => ({
  categoricalSelectionState: state.categoricalSelectionState
}))
class Categories extends React.Component {
  render() {
    const { categoricalSelectionState } = this.props;
    if (!categoricalSelectionState) return null;

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
        {_.map(categoricalSelectionState, (catState, catName) => (
          <Category key={catName} metadataField={catName} />
        ))}
      </div>
    );
  }
}

export default Categories;
