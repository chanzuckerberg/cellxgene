import React from "react";
import _ from "lodash";
import { connect } from "react-redux";
import { Flipper, Flipped } from "react-flip-toolkit";

import * as globals from "../../../globals";
import Value from "../value";

@connect((state) => ({
  categoricalSelection: state.categoricalSelection,
}))
class Category extends React.Component {
  constructor(props) {
    super(props);
    this.state = {};
  }

  renderCategoryItems(optTuples) {
    const { metadataField, isUserAnno } = this.props;

    return _.map(optTuples, (tuple, i) => {
      return (
        <Flipped key={tuple[1]} flipId={tuple[1]}>
          {(flippedProps) => (
            <Value
              isUserAnno={isUserAnno}
              optTuples={optTuples}
              key={tuple[1]}
              metadataField={metadataField}
              categoryIndex={tuple[1]}
              i={i}
              flippedProps={flippedProps}
            />
          )}
        </Flipped>
      );
    });
  }

  render() {
    const {
      metadataField,
      categoricalSelection,
      children,
      isExpanded,
    } = this.props;
    const { isTruncated } = categoricalSelection[metadataField];
    const cat = categoricalSelection[metadataField];
    const optTuples = [...cat.categoryValueIndices];
    const optTuplesAsKey = _.map(optTuples, (t) => t[0]).join(""); // animation

    return (
      <div
        style={{
          maxWidth: globals.maxControlsWidth,
        }}
        data-testclass="category"
        data-testid={`category-${metadataField}`}
      >
        <div
          style={{
            display: "flex",
            justifyContent: "space-between",
            alignItems: "baseline",
          }}
        >
          {children}
        </div>
        <div style={{ marginLeft: 26 }}>
          <Flipper spring="veryGentle" flipKey={optTuplesAsKey}>
            {isExpanded ? this.renderCategoryItems(optTuples) : null}
          </Flipper>
        </div>
        <div>
          {isExpanded && isTruncated ? (
            <p style={{ paddingLeft: 15 }}>... truncated list ...</p>
          ) : null}
        </div>
      </div>
    );
  }
}

export default Category;
