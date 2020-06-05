import React from "react";
import { Flipper, Flipped } from "react-flip-toolkit";

import * as globals from "../../../globals";
import Value from "../value";

class Category extends React.Component {
  constructor(props) {
    super(props);
    this.state = {};
  }

  renderCategoryItems(optTuples) {
    const { metadataField, isUserAnno, categorySummary } = this.props;

    return optTuples.map((tuple, i) => {
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
              categorySummary={categorySummary}
            />
          )}
        </Flipped>
      );
    });
  }

  render() {
    const { metadataField, categorySummary, children, isExpanded } = this.props;
    const { isTruncated } = categorySummary;
    const optTuples = [...categorySummary.categoryValueIndices];
    const optTuplesAsKey = optTuples.map((t) => t[0]).join(""); // animation

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
