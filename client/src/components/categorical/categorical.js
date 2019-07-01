// jshint esversion: 6
import React from "react";
import _ from "lodash";
import { connect } from "react-redux";
import * as globals from "../../globals";
import Category from "./category";

@connect(state => ({
  categoricalSelection: state.categoricalSelection
}))
class Categories extends React.Component {
  render() {
    const { categoricalSelection } = this.props;
    if (!categoricalSelection) return null;

    return (
      <div
        style={{
          padding: globals.leftSidebarSectionPadding
        }}
      >
        {_.map(categoricalSelection, (catState, catName) => (
          <Category key={catName} metadataField={catName} />
        ))}
      </div>
    );
  }
}

export default Categories;
